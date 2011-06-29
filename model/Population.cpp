/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include "inputData.h"
#include "Monitoring/Surveys.h"
#include "Monitoring/Continuous.h"

#include "Transmission/TransmissionModel.h"
#include "Transmission/Vector/VectorTransmission.h"

#include "Host/Human.h"
#include "Host/NeonatalMortality.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/CaseManagementCommon.h"
#include "Pathogenesis/PathogenesisModel.h"
#include "PkPd/PkPdModel.h"

#include "util/errors.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"

#include <cmath>
#include <boost/format.hpp>
#include <boost/assign.hpp>

namespace OM
{
    using namespace OM::util;
    using namespace boost::assign;

// -----  Population: static data / methods  -----

void Population::init()
{
    Host::Human::initHumanParameters();
    Host::NeonatalMortality::init();
    PkPd::PkPdModel::init();
    
    AgeStructure::init ();
}

void Population::clear()
{
    PkPd::PkPdModel::cleanup ();
    Host::Human::clear();
}

void Population::staticCheckpoint (istream& stream)
{
    Host::NeonatalMortality::staticCheckpoint (stream);
    Clinical::ClinicalModel::staticCheckpoint (stream);
    PkPd::PkPdModel::staticCheckpoint (stream);
}
void Population::staticCheckpoint (ostream& stream)
{
    Host::NeonatalMortality::staticCheckpoint (stream);
    Clinical::ClinicalModel::staticCheckpoint (stream);
    PkPd::PkPdModel::staticCheckpoint (stream);
}


// -----  non-static methods: creation/destruction, checkpointing  -----

Population::Population()
    : populationSize (InputData().getDemography().getPopSize())
{
    using Monitoring::Continuous;
    Continuous::registerCallback( "hosts", "\thosts", MakeDelegate( this, &Population::ctsHosts ) );
    // Age groups are currently hard-coded.
    ctsDemogAgeGroups += 1.0, 5.0, 10.0, 15.0, 25.0;
    ostringstream ctsDemogTitle;
    BOOST_FOREACH( double ubound, ctsDemogAgeGroups ){
        ctsDemogTitle << "\thost % ≤ " << ubound;
    }
    Continuous::registerCallback( "host demography", ctsDemogTitle.str(), MakeDelegate( this, &Population::ctsHostDemography ) );
    Continuous::registerCallback( "recent births", "\trecent births", MakeDelegate( this, &Population::ctsRecentBirths ) );
    Continuous::registerCallback( "patent hosts", "\tpatent hosts", MakeDelegate( this, &Population::ctsPatentHosts ) );
    Continuous::registerCallback( "immunity h", "\timmunity h", MakeDelegate( this, &Population::ctsImmunityh ) );
    Continuous::registerCallback( "immunity Y", "\timmunity Y", MakeDelegate( this, &Population::ctsImmunityY ) );
    Continuous::registerCallback( "median immunity Y", "\tmedian immunity Y", MakeDelegate( this, &Population::ctsMedianImmunityY ) );
    Continuous::registerCallback( "human age availability", "\thuman age availability", MakeDelegate( this, &Population::ctsMeanAgeAvailEffect ) );
    Continuous::registerCallback( "nets owned", "\tnets owned", MakeDelegate( this, &Population::ctsNetsOwned ) );
    Continuous::registerCallback( "mean hole index", "\tmean hole index", MakeDelegate( this, &Population::ctsNetHoleIndex ) );
    Continuous::registerCallback( "mean insecticide content", "\tmean insecticide content", MakeDelegate( this, &Population::ctsNetInsecticideContent ) );
    
    _transmissionModel = Transmission::TransmissionModel::createTransmissionModel(populationSize);
}

Population::~Population()
{
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        iter->destroy();
    }
    delete _transmissionModel;
}

void Population::checkpoint (istream& stream)
{
    size_t popSize; // must be type of population.size()
    popSize & stream;
    if (popSize > size_t (populationSize))
        throw util::checkpoint_error( (boost::format("pop size (%1%) exceeds that given in scenario.xml") %popSize).str() );
    for (size_t i = 0; i < popSize && !stream.eof(); ++i) {
        // Note: calling this constructor of Host::Human is slightly wasteful, but avoids the need for another
        // ctor and leaves less opportunity for uninitialized memory.
        population.push_back (Host::Human (*_transmissionModel, TimeStep(0)));
        population.back() & stream;
    }
    if (population.size() != popSize)
        throw util::checkpoint_error( (boost::format("Population: out of data (read %1% humans)") %population.size() ).str() );
}
void Population::checkpoint (ostream& stream)
{
    population.size() & stream;
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter)
        (*iter) & stream;
}

void Population::preMainSimInit ()
{
    Host::InfectionIncidenceModel::initMainSimulation();
    Clinical::ClinicalModel::initMainSimulation();
    recentBirths = 0;
}

void Population::createInitialHumans ()
{
    /* We create a whole population here, regardless of whether humans can
    survive until start of vector init (vector model needs a whole population
    structure in any case). However, we don't update humans known not to survive
    until vector init, which saves computation and memory (no infections). */
    
    int cumulativePop = 0;
    for (TimeStep iage = AgeStructure::getMaxTimestepsPerLife() - TimeStep(1);
         iage >= TimeStep(0); --iage )
    {
	int targetPop = AgeStructure::targetCumPop (iage, populationSize);
	while (cumulativePop < targetPop) {
	    newHuman (-iage);
	    ++cumulativePop;
	}
    }
    
    // Vector setup dependant on human population structure (we *want* to
    // include all humans, whether they'll survive to vector init phase or not).
    _transmissionModel->setupNv0 (population, populationSize);
}


// -----  non-static methods: simulation loop  -----

void Population::newHuman (TimeStep dob)
{
    util::streamValidate( dob.asInt() );
    population.push_back (Host::Human (*_transmissionModel, dob));
    ++recentBirths;
}

void Population::update1()
{
    // This should only use humans being updated: otherwise too small a proportion
    // will be infected. However, we don't have another number to use instead.
    // TODO: add some value to be used until old-enough humans are updated?
    Host::NeonatalMortality::update (population);
    
    // Should be called before Human::update().
    // Vector: This needs the whole population (it is an approximation before all humans are updated).
    // NonVector: Doesn't matter whether non-updated humans are included (value isn't used
    // before all humans are updated).
    _transmissionModel->update (population, populationSize);

    //NOTE: other parts of code are not set up to handle changing population size. Also
    // populationSize is assumed to be the _actual and exact_ population size by other code.
    //targetPop is the population size at time t allowing population growth
    //int targetPop = (int) (populationSize * exp (AgeStructure::rho * TimeStep::simulation));
    int targetPop = populationSize;
    int cumPop = 0;

    // Update each human in turn
    //std::cout<<" time " <<t<<std::endl;
    HumanIter last = population.end();
    --last;
    for (HumanIter iter = population.begin(); iter != population.end();) {
        // Update human, and remove if too old:
        if (iter->update (*this, _transmissionModel,
		/* Only include humans who can survive until vector init.
		Note: we could exclude more humans due to age distribution,
		but how many extra to leave due to deaths isn't obvious. */
		TimeStep::intervalsPerYear + iter->getDateOfBirth() > TimeStep(0)
	    )) {
            iter->destroy();
            iter = population.erase (iter);
            continue;
        }

        //BEGIN Population size & age structure
        ++cumPop;
        TimeStep age = (TimeStep::simulation - iter->getDateOfBirth());

        // if (Actual number of people so far > target population size for this age) ...
        if (cumPop > AgeStructure::targetCumPop (age, targetPop)) {
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
        newHuman (TimeStep::simulation);
        //++nCounter;
        ++cumPop;
    }
    
    _transmissionModel->updateSummaries();
}


// -----  non-static methods: reporting  -----

void Population::ctsHosts (ostream& stream){
    // this option is intended for debugging human initialization; normally this should equal populationSize.
    stream << '\t' << population.size();
}
void Population::ctsHostDemography (ostream& stream){
    list<Host::Human>::reverse_iterator it = population.rbegin();
    int cumCount = 0;
    BOOST_FOREACH( double ubound, ctsDemogAgeGroups ){
        while( it != population.rend() && it->getAgeInYears() < ubound ){
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
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        if( iter->getWithinHostModel().parasiteDensityDetectible() )
            ++patent;
    }
    stream << '\t' << patent;
}
void Population::ctsImmunityh (ostream& stream){
    double x = 0.0;
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        x += iter->getWithinHostModel().getCumulativeh();
    }
    x /= populationSize;
    stream << '\t' << x;
}
void Population::ctsImmunityY (ostream& stream){
    double x = 0.0;
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        x += iter->getWithinHostModel().getCumulativeY();
    }
    x /= populationSize;
    stream << '\t' << x;
}
void Population::ctsMedianImmunityY (ostream& stream){
    vector<double> list;
    list.reserve( populationSize );
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        list.push_back( iter->getWithinHostModel().getCumulativeY() );
    }
    sort( list.begin(), list.end() );
    double x;
    if( populationSize % 2 == 0 ){
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
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        if( !iter->perHostTransmission.isOutsideTransmission() ){
            ++nHumans;
            avail += iter->perHostTransmission.relativeAvailabilityAge(iter->getAgeInYears());
        }
    }
    stream << '\t' << avail/nHumans;
}
void Population::ctsNetsOwned (ostream& stream){
    int nNets = 0;
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        if( iter->perHostTransmission.getITN().timeOfDeployment() >= TimeStep(0) )
            ++nNets;
    }
    stream << '\t' << nNets;
}
void Population::ctsNetHoleIndex (ostream& stream){
    double meanVar = 0.0;
    int nNets = 0;
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        if( iter->perHostTransmission.getITN().timeOfDeployment() >= TimeStep(0) ){
            ++nNets;
            meanVar += iter->perHostTransmission.getITN().getHoleIndex();
        }
    }
    stream << '\t' << meanVar/nNets;
}
void Population::ctsNetInsecticideContent (ostream& stream){
    const Transmission::VectorTransmission* vt = reinterpret_cast<Transmission::VectorTransmission*>(_transmissionModel);
    if( vt == 0 ){
        throw util::xml_scenario_error("mean insecticide content: invalid in non-vector mode!");
    }
    const Transmission::ITNParams& params = vt->getITNParams();
    double meanVar = 0.0;
    int nNets = 0;
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        if( iter->perHostTransmission.getITN().timeOfDeployment() >= TimeStep(0) ){
            ++nNets;
            meanVar += iter->perHostTransmission.getITN().getInsecticideContent(params);
        }
    }
    stream << '\t' << meanVar/nNets;
}

void Population::newSurvey ()
{
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        iter->summarize();
    }
    _transmissionModel->summarize( *Monitoring::Surveys.current );
}

void Population::flushReports (){
    for (HumanIter iter = population.begin(); iter != population.end(); ++iter) {
        iter->flushReports();
    }
}    

}
