/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "Population.h"
#include "mon/Continuous.h"

#include "Host/Human.h"
#include "Host/NeonatalMortality.h"
#include "Host/WithinHost/WHInterface.h"
#include "Host/WithinHost/Genotypes.h"
#include "Host/WithinHost/Diagnostic.h"
#include "Clinical/ClinicalModel.h"
#include "Transmission/TransmissionModel.h"

#include "util/errors.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include <schema/scenario.h>

#include <cmath>

namespace OM
{
    using namespace std::placeholders;
    using namespace OM::util;
    using Transmission::TransmissionModel;

// -----  non-static methods: creation/destruction, checkpointing  -----

Population::Population(size_t populationSize) : populationSize (populationSize){}

namespace population
{
    Population *createPopulation(size_t populationSize)
    {
        Population *population = new Population(populationSize);

        ostringstream ctsDemogTitle;
        for( double ubound : population->ctsDemogAgeGroups )
            ctsDemogTitle << "\thost % ≤ " << ubound;

        mon::Continuous.registerCallback( "hosts", "\thosts", &population::ctsHosts );
        mon::Continuous.registerCallback( "host demography", ctsDemogTitle.str(), &population::ctsHostDemography);
        mon::Continuous.registerCallback( "recent births", "\trecent births", &population::ctsRecentBirths);
        mon::Continuous.registerCallback( "patent hosts", "\tpatent hosts", &population::ctsPatentHosts);
        mon::Continuous.registerCallback( "immunity h", "\timmunity h", &population::ctsImmunityh);
        mon::Continuous.registerCallback( "immunity Y", "\timmunity Y", &population::ctsImmunityY);
        mon::Continuous.registerCallback( "median immunity Y", "\tmedian immunity Y", &population::ctsMedianImmunityY);
        mon::Continuous.registerCallback( "human age availability", "\thuman age availability", &population::ctsMeanAgeAvailEffect);
        mon::Continuous.registerCallback( "ITN coverage", "\tITN coverage", &population::ctsITNCoverage);
        mon::Continuous.registerCallback( "IRS coverage", "\tIRS coverage", &population::ctsIRSCoverage);
        mon::Continuous.registerCallback( "GVI coverage", "\tGVI coverage", &population::ctsGVICoverage);

        return population;
    }

    void checkpoint (Population &population, istream& stream)
    {
        population.populationSize & stream;
        population.recentBirths & stream;

        auto &humans = population.humans;
        
        for(size_t i = 0; i < population.populationSize && !stream.eof(); ++i) {
            // Note: calling this constructor of Host::Human is slightly wasteful, but avoids the need for another
            // ctor and leaves less opportunity for uninitialized memory.
            humans.push_back( Host::Human (sim::zero()) );
            humans.back() & stream;
        }
        if (humans.size() != population.populationSize)
            throw util::checkpoint_error("Population: out of data (read " + to_string(humans.size()) + " humans)");
    }
    void checkpoint (Population &population, ostream& stream)
    {
        population.populationSize & stream;
        population.recentBirths & stream;
        auto &humans = population.humans;

        for(Host::Human &human : humans)
            human & stream;
    }

    void createInitialHumans(Population &population)
    {
        const int populationSize = population.populationSize;
        auto &humans = population.humans;

        /* We create a whole population here, regardless of whether humans can
        survive until start of vector init (vector model needs a whole population
        structure in any case). However, we don't update humans known not to survive
        until vector init, which saves computation and memory (no infections). */
        
        int cumulativePop = 0;
        for(size_t iage_prev = AgeStructure::getMaxTStepsPerLife(), iage = iage_prev - 1;
             iage_prev > 0; iage_prev = iage, iage -= 1 )
        {
            int targetPop = AgeStructure::targetCumPop( iage, populationSize );
            while (cumulativePop < targetPop) {
                SimTime dob = sim::zero() - sim::fromTS(iage);
                util::streamValidate( dob );
                humans.push_back( Host::Human (dob) );
                ++cumulativePop;
            }
        }
        
        // Vector setup dependant on human population structure (we *want* to
        // include all humans, whether they'll survive to vector init phase or not).
        assert( sim::now() == sim::zero() );      // assumed below
    }


    // -----  non-static methods: simulation loop  -----

    void update(Population &population, Transmission::TransmissionModel& transmission, SimTime firstVecInitTS ){
        const int populationSize = population.populationSize;
        int &recentBirths = population.recentBirths;
        auto &humans = population.humans;

        // This should only use humans being updated: otherwise too small a proportion
        // will be infected. However, we don't have another number to use instead.
        // NOTE: no neonatal mortalities will occur in the first 20 years of warmup
        // (until humans old enough to be pregnate get updated and can be infected).
        Host::NeonatalMortality::update (humans);
        
        // Update each human in turn
        for (Host::Human& human : humans) {
            // Update human, and remove if too old.
            // We only need to update humans who will survive past the end of the
            // "one life span" init phase (this is an optimisation). lastPossibleTS
            // is the time step they die at (some code still runs on this step).
            SimTime lastPossibleTS = human.getDateOfBirth() + sim::maxHumanAge();   // this is last time of possible update
            if (lastPossibleTS >= firstVecInitTS)
                human.update(transmission);
        }
        
        //NOTE: other parts of code are not set up to handle changing population size. Also
        // populationSize is assumed to be the _actual and exact_ population size by other code.
        //targetPop is the population size at time t allowing population growth
        //int targetPop = (int) (populationSize * exp( AgeStructure::rho * sim::ts1().inSteps() ));
        int targetPop = populationSize;
        int cumPop = 0;

        for (auto it = humans.begin(); it != humans.end();) {
            bool isDead = it->remove();
            // if (Actual number of people so far > target population size for this age)
            // "outmigrate" some to maintain population shape
            //NOTE: better to use age(sim::ts0())? Possibly, but the difference will not be very significant.
            // Also see targetPop = ... comment above
            bool outmigrate = cumPop >= AgeStructure::targetCumPop(sim::inSteps(it->age(sim::ts1())), targetPop);
            
            if( isDead || outmigrate ){
                it = humans.erase (it);
                continue;
            }
            ++cumPop;
            ++it;
        } // end of per-human updates

        // increase population size to targetPop
        recentBirths += (targetPop - cumPop);
        while (cumPop < targetPop) {
            // humans born at end of this time step = beginning of next, hence ts1
            humans.push_back( Host::Human (sim::ts1()) );
            ++cumPop;
        }
    }


    // -----  non-static methods: reporting  -----

    void newSurvey (Population &population)
    {
        for(Host::Human &human : population.humans) {
            human.summarize();
        }
    }

    void flushReports (Population &population)
    {
        for(Host::Human &human : population.humans) {
            human.flushReports();
        }
    }  

    void ctsHosts (Population &population, ostream& stream){
        // this option is intended for debugging human initialization; normally this should equal populationSize.
        stream << '\t' << population.populationSize;
    }

    void ctsHostDemography (Population &population, ostream& stream){
        auto iter = population.humans.crbegin();
        int cumCount = 0;
        for( double ubound : population.ctsDemogAgeGroups ){
            while( iter != population.humans.crend() && sim::inYears(iter->age(sim::now())) < ubound ){
                ++cumCount;
                ++iter;
            }
            stream << '\t' << cumCount;
        }
    }

    void ctsRecentBirths (Population &population, ostream& stream){
        stream << '\t' << population.recentBirths;
        population.recentBirths = 0;
    }

    void ctsPatentHosts (Population &population, ostream& stream){
        int patent = 0;
        for(Host::Human &human : population.humans) {
            auto diag = WithinHost::diagnostics::monitoringDiagnostic();
            if( human.getWithinHostModel().diagnosticResult(human.rng(), diag) )
                ++patent;
        }
        stream << '\t' << patent;
    }

    void ctsImmunityh (Population &population, ostream& stream){
        double x = 0.0;
        for(const Host::Human &human : population.humans) {
            x += human.getWithinHostModel().getCumulative_h();
        }
        x /= population.populationSize;
        stream << '\t' << x;
    }

    void ctsImmunityY (Population &population, ostream& stream){
        double x = 0.0;
        for(const Host::Human &human : population.humans) {
            x += human.getWithinHostModel().getCumulative_Y();
        }
        x /= population.populationSize;
        stream << '\t' << x;
    }

    void ctsMedianImmunityY (Population &population, ostream& stream){
        vector<double> list;
        list.reserve( population.populationSize );
        for(const Host::Human &human : population.humans) {
            list.push_back( human.getWithinHostModel().getCumulative_Y() );
        }
        sort( list.begin(), list.end() );
        double x;
        if( mod_nn(population.populationSize, 2) == 0 ){
            size_t i = population.populationSize / 2;
            x = (list[i-1]+list[i])/2.0;
        }else{
            x = list[population.populationSize / 2];
        }
        stream << '\t' << x;
    }

    void ctsMeanAgeAvailEffect (Population &population, ostream& stream){
        int nHumans = 0;
        double avail = 0.0;
        for(Host::Human &human : population.humans) {
            if( !human.perHostTransmission.isOutsideTransmission() ){
                ++nHumans;
                avail += human.perHostTransmission.relativeAvailabilityAge(sim::inYears(human.age(sim::now())));
            }
        }
        stream << '\t' << avail/nHumans;
    }

    void ctsITNCoverage (Population &population, ostream& stream){
        int nActive = 0;
        for(const Host::Human &human : population.humans) {
            nActive += human.perHostTransmission.hasActiveInterv( interventions::Component::ITN );
        }
        double coverage = static_cast<double>(nActive) / population.populationSize;
        stream << '\t' << coverage;
    }

    void ctsIRSCoverage (Population &population, ostream& stream){
        int nActive = 0;
        for(const Host::Human &human : population.humans) {
            nActive += human.perHostTransmission.hasActiveInterv( interventions::Component::IRS );
        }
        double coverage = static_cast<double>(nActive) / population.populationSize;
        stream << '\t' << coverage;
    }
    
    void ctsGVICoverage (Population &population, ostream& stream){
        int nActive = 0;
        for(const Host::Human &human : population.humans) {
            nActive += human.perHostTransmission.hasActiveInterv( interventions::Component::GVI );
        }
        double coverage = static_cast<double>(nActive) / population.populationSize;
        stream << '\t' << coverage;
    }
}

}

