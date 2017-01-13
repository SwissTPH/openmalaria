/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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
#include "Monitoring/Continuous.h"

#include "Host/Human.h"
#include "Host/NeonatalMortality.h"
#include "WithinHost/WHInterface.h"
#include "WithinHost/Genotypes.h"
#include "WithinHost/Diagnostic.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/CaseManagementCommon.h"
#include "Transmission/TransmissionModel.h"

#include "util/errors.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include <schema/scenario.h>

#include <cmath>
#include <boost/format.hpp>
#include <boost/assign.hpp>

namespace OM
{
    using namespace OM::util;
    using namespace boost::assign;
    using Transmission::TransmissionModel;

// -----  Population: static data / methods  -----

void Population::init( const Parameters& parameters, const scnXml::Scenario& scenario )
{
    Host::Human::init( parameters, scenario );
    Host::NeonatalMortality::init( scenario.getModel().getClinical() );
    
    AgeStructure::init( scenario.getDemography() );
}

void Population::staticCheckpoint (istream& stream)
{
    Host::NeonatalMortality::staticCheckpoint (stream);
    Host::InfectionIncidenceModel::staticCheckpoint (stream);
    Clinical::staticCheckpointCMCommon (stream);
    WithinHost::Genotypes::staticCheckpoint(stream);
}
void Population::staticCheckpoint (ostream& stream)
{
    Host::NeonatalMortality::staticCheckpoint (stream);
    Host::InfectionIncidenceModel::staticCheckpoint (stream);
    Clinical::staticCheckpointCMCommon (stream);
    WithinHost::Genotypes::staticCheckpoint(stream);
}


// -----  non-static methods: creation/destruction, checkpointing  -----

Population::Population(size_t populationSize)
    : populationSize (populationSize), recentBirths(0)
{
    using Monitoring::Continuous;
    Continuous.registerCallback( "hosts", "\thosts", MakeDelegate( this, &Population::ctsHosts ) );
    // Age groups are currently hard-coded.
    ctsDemogAgeGroups += 1.0, 5.0, 10.0, 15.0, 25.0;
    ostringstream ctsDemogTitle;
    foreach( double ubound, ctsDemogAgeGroups ){
        ctsDemogTitle << "\thost % ≤ " << ubound;
    }
    Continuous.registerCallback( "host demography", ctsDemogTitle.str(),
        MakeDelegate( this, &Population::ctsHostDemography ) );
    Continuous.registerCallback( "recent births", "\trecent births",
        MakeDelegate( this, &Population::ctsRecentBirths ) );
    Continuous.registerCallback( "patent hosts", "\tpatent hosts",
        MakeDelegate( this, &Population::ctsPatentHosts ) );
    Continuous.registerCallback( "immunity h", "\timmunity h",
        MakeDelegate( this, &Population::ctsImmunityh ) );
    Continuous.registerCallback( "immunity Y", "\timmunity Y",
        MakeDelegate( this, &Population::ctsImmunityY ) );
    Continuous.registerCallback( "median immunity Y", "\tmedian immunity Y",
        MakeDelegate( this, &Population::ctsMedianImmunityY ) );
    Continuous.registerCallback( "human age availability",
        "\thuman age availability",
        MakeDelegate( this, &Population::ctsMeanAgeAvailEffect ) );
    Continuous.registerCallback( "ITN coverage", "\tITN coverage",
        MakeDelegate( this, &Population::ctsITNCoverage ) );
    Continuous.registerCallback( "IRS coverage", "\tIRS coverage",
        MakeDelegate( this, &Population::ctsIRSCoverage ) );
    Continuous.registerCallback( "GVI coverage", "\tGVI coverage",
        MakeDelegate( this, &Population::ctsGVICoverage ) );
    // "nets owned" replaced by "ITN coverage"
//     Continuous.registerCallback( "nets owned", "\tnets owned",
//         MakeDelegate( this, &Population::ctsNetsOwned ) );
//     Continuous.registerCallback( "mean hole index", "\tmean hole index",
//         MakeDelegate( this, &Population::ctsNetHoleIndex ) );
}

Population::~Population()
{
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        iter->destroy();
    }
}

void Population::checkpoint (istream& stream)
{
    populationSize & stream;
    recentBirths & stream;
    
    for(size_t i = 0; i < populationSize && !stream.eof(); ++i) {
        // Note: calling this constructor of Host::Human is slightly wasteful, but avoids the need for another
        // ctor and leaves less opportunity for uninitialized memory.
        population.push_back( new Host::Human (SimTime::zero()) );
        population.back() & stream;
    }
    if (population.size() != populationSize)
        throw util::checkpoint_error(
            (boost::format("Population: out of data (read %1% humans)") %population.size() ).str() );
}
void Population::checkpoint (ostream& stream)
{
    populationSize & stream;
    recentBirths & stream;
    
    for(Iter iter = population.begin(); iter != population.end(); ++iter)
        (*iter) & stream;
}

void Population::preMainSimInit ()
{
    Host::InfectionIncidenceModel::initMainSimulation();
    Clinical::mainSimInitCMCommon();
    WithinHost::Genotypes::startMainSim();
    recentBirths = 0;
}

void Population::createInitialHumans()
{
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
            newHuman( SimTime::zero() - SimTime::fromTS(iage) );
            ++cumulativePop;
        }
    }
    
    // Vector setup dependant on human population structure (we *want* to
    // include all humans, whether they'll survive to vector init phase or not).
    assert( sim::now() == SimTime::zero() );      // assumed below
}


// -----  non-static methods: simulation loop  -----

void Population::newHuman( SimTime dob ){
    util::streamValidate( dob.raw() );
    population.push_back( new Host::Human (dob) );
    ++recentBirths;
}

void Population::update1( SimTime firstVecInitTS ){
    // This should only use humans being updated: otherwise too small a proportion
    // will be infected. However, we don't have another number to use instead.
    // NOTE: no neonatal mortalities will occur in the first 20 years of warmup
    // (until humans old enough to be pregnate get updated and can be infected).
    Host::NeonatalMortality::update (*this);
    
    //NOTE: other parts of code are not set up to handle changing population size. Also
    // populationSize is assumed to be the _actual and exact_ population size by other code.
    //targetPop is the population size at time t allowing population growth
    //int targetPop = (int) (populationSize * exp( AgeStructure::rho * sim::ts1().inSteps() ));
    int targetPop = populationSize;
    int cumPop = 0;

    // Update each human in turn
    //std::cout<<" time " <<t<<std::endl;
    Iter last = population.end();
    --last;
    for(Iter iter = population.begin(); iter != population.end();) {
        // Update human, and remove if too old.
        // We only need to update humans who will survive past the end of the
        // "one life span" init phase (this is an optimisation). lastPossibleTS
        // is the time step they die at (some code still runs on this step).
        SimTime lastPossibleTS = iter->getDateOfBirth() + sim::maxHumanAge();   // this is last time of possible update
        bool updateHuman = lastPossibleTS >= firstVecInitTS;
        bool isDead = iter->update(updateHuman);
        if( isDead ){
            iter->destroy();
            iter = population.erase (iter);
            continue;
        }
        
        //BEGIN Population size & age structure
        ++cumPop;

        // if (Actual number of people so far > target population size for this age)
        // "outmigrate" some to maintain population shape
        //NOTE: better to use age(sim::ts0())? Possibly, but the difference will not be very significant.
        // Also see targetPop = ... comment above
        if( cumPop > AgeStructure::targetCumPop(iter->age(sim::ts1()).inSteps(), targetPop) ){
            --cumPop;
            iter->destroy();
            iter = population.erase (iter);
            continue;
        }
        //END Population size & age structure
        ++iter;
    } // end of per-human updates

    // increase population size to targetPop
    while (cumPop < targetPop) {
        // humans born at end of this time step = beginning of next, hence ts1
        newHuman( sim::ts1() );
        //++nCounter;
        ++cumPop;
    }
}


// -----  non-static methods: reporting  -----

void Population::ctsHosts (ostream& stream){
    // this option is intended for debugging human initialization; normally this should equal populationSize.
    stream << '\t' << population.size();
}
void Population::ctsHostDemography (ostream& stream){
    Population::ConstReverseIter it = population.crbegin();
    int cumCount = 0;
    foreach( double ubound, ctsDemogAgeGroups ){
        while( it != population.crend() && it->age(sim::now()).inYears() < ubound ){
            ++cumCount;
            ++it;
        }
        stream << '\t' << cumCount;
    }
}
void Population::ctsRecentBirths (ostream& stream){
    stream << '\t' << recentBirths;
    recentBirths = 0;
}
void Population::ctsPatentHosts (ostream& stream){
    int patent = 0;
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        if( iter->getWithinHostModel().diagnosticResult(WithinHost::diagnostics::monitoringDiagnostic()) )
            ++patent;
    }
    stream << '\t' << patent;
}
void Population::ctsImmunityh (ostream& stream){
    double x = 0.0;
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        x += iter->getWithinHostModel().getCumulative_h();
    }
    x /= populationSize;
    stream << '\t' << x;
}
void Population::ctsImmunityY (ostream& stream){
    double x = 0.0;
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        x += iter->getWithinHostModel().getCumulative_Y();
    }
    x /= populationSize;
    stream << '\t' << x;
}
void Population::ctsMedianImmunityY (ostream& stream){
    vector<double> list;
    list.reserve( populationSize );
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        list.push_back( iter->getWithinHostModel().getCumulative_Y() );
    }
    sort( list.begin(), list.end() );
    double x;
    if( mod_nn(populationSize, 2) == 0 ){
        size_t i = populationSize / 2;
        x = (list[i-1]+list[i])/2.0;
    }else{
        x = list[populationSize / 2];
    }
    stream << '\t' << x;
}
void Population::ctsMeanAgeAvailEffect (ostream& stream){
    int nHumans = 0;
    double avail = 0.0;
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        if( !iter->perHostTransmission.isOutsideTransmission() ){
            ++nHumans;
            avail += iter->perHostTransmission.relativeAvailabilityAge(iter->age(sim::now()).inYears());
        }
    }
    stream << '\t' << avail/nHumans;
}
void Population::ctsITNCoverage (ostream& stream){
    int nActive = 0;
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        nActive += iter->perHostTransmission.hasActiveInterv( interventions::Component::ITN );
    }
    double coverage = static_cast<double>(nActive) / populationSize;
    stream << '\t' << coverage;
}
void Population::ctsIRSCoverage (ostream& stream){
    int nActive = 0;
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        nActive += iter->perHostTransmission.hasActiveInterv( interventions::Component::IRS );
    }
    double coverage = static_cast<double>(nActive) / populationSize;
    stream << '\t' << coverage;
}
void Population::ctsGVICoverage (ostream& stream){
    int nActive = 0;
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        nActive += iter->perHostTransmission.hasActiveInterv( interventions::Component::GVI );
    }
    double coverage = static_cast<double>(nActive) / populationSize;
    stream << '\t' << coverage;
}
// void Population::ctsNetHoleIndex (ostream& stream){
//     double meanVar = 0.0;
//     int nNets = 0;
//     for(Iter iter = population.begin(); iter != population.end(); ++iter) {
//         if( iter->perHostTransmission.getITN().timeOfDeployment() >= SimTime::zero() ){
//             ++nNets;
//             meanVar += iter->perHostTransmission.getITN().getHoleIndex();
//         }
//     }
//     stream << '\t' << meanVar/nNets;
// }

void Population::newSurvey ()
{
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        iter->summarize();
    }
}

void Population::flushReports (){
    for(Iter iter = population.begin(); iter != population.end(); ++iter) {
        iter->flushReports();
    }
}    

}

