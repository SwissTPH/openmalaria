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

#include "Transmission/NonVectorModel.h"
#include "Transmission/PerHost.h"
#include "Host/Human.h"
#include "WithinHost/Genotypes.h"
#include "mon/info.h"
#include "util/random.h"
#include "util/vectors.h"
#include "util/StreamValidator.h"
#include "util/checkpoint_containers.h"
#include <limits>
#include <cmath>

namespace OM { namespace Transmission {
    namespace vectors = util::vectors;

//static (class) variables
const double NonVectorModel::totalInfectionrateVariance= 1.0;
const double NonVectorModel::min_EIR_mult= 0.01;

const int nYearsWarmupData = 5;

NonVectorModel::NonVectorModel(const scnXml::Entomology& entoData,
        const scnXml::NonVector& nonVectorData) :
    TransmissionModel(entoData, 1/*this model doesn't support multiple genotypes*/),
    nSpore( sim::fromDays( nonVectorData.getEipDuration() ) )
{
    laggedKappa.resize( nSpore.inSteps() + 1, 0.0 );
    
    vector<int> nDays (sim::stepsPerYear(), 0);
    //The minimum EIR allowed in the array. The product of the average EIR and a constant.
    double minEIR=min_EIR_mult*averageEIR(nonVectorData);
    
    const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();
    if( daily.size() < static_cast<size_t>(sim::oneYear().inDays()) )
        throw util::xml_scenario_error( "insufficient EIRDaily data for a year" );
    
    for( SimTime mpcday = sim::zero(), endDay = sim::fromDays(daily.size());
         mpcday < endDay; mpcday += sim::oneDay() )
    {
        double EIRdaily = std::max(static_cast<double>(daily[mpcday.inDays()]), minEIR);
        
        // Index 0 of initialisationEIR refers to the EIR affecting the
        // first day(s) of the year. Correspondingly, the first 1 or 5 values
        // of EIRDaily affect this (1- or 5-day) time-step.
        size_t i = mod_nn(mpcday.inSteps(), sim::stepsPerYear());
        
        nDays[i] += 1;
        initialisationEIR[i] += EIRdaily;
    }
    
    // Calculate total annual EIR
    // divide by number of records assigned to each interval (usually one per day)
    for( size_t indTS = 0; indTS < sim::stepsPerYear(); indTS += 1 ){
        initialisationEIR[indTS] *= sim::oneTS().inDays() / static_cast<double>( nDays[indTS] );
        annualEIR += initialisationEIR[indTS];
    }
    
    initialKappa.assign( sim::fromYearsI(nYearsWarmupData).inSteps(), 0.0 );
}

NonVectorModel::~NonVectorModel () {}

void NonVectorModel::init2 (const Population&) {
    // no set-up needed; just indicate we're ready to roll:
    simulationMode = forcedEIR;
}

const char* viError = "vector model interventions can not be used with the non-vector model";
void NonVectorModel::initVectorInterv( const scnXml::Description::AnophelesSequence& list,
        size_t instance, const string& name )
{
    throw util::xml_scenario_error( viError );
}
void NonVectorModel::initVectorTrap( const scnXml::VectorTrap::DescriptionSequence list,
        size_t instance, const scnXml::VectorTrap::NameOptional name )
{
    throw util::xml_scenario_error( viError );
}


void NonVectorModel::scaleEIR (double factor){
    vectors::scale( initialisationEIR, factor );
    annualEIR = vectors::sum( initialisationEIR );
}
#if 0
void NonVectorModel::scaleXML_EIR (scnXml::Entomology& ed, double factor) const{
    assert( ed.getNonVector().present() );
    scnXml::NonVector::EIRDailySequence& daily = ed.getNonVector().get().getEIRDaily();
    
    for( scnXml::NonVector::EIRDailyIterator it = daily.begin();
	it != daily.end(); ++it ){
	double old = *it;
	*it = old * factor;
    }
}
#endif


SimTime NonVectorModel::minPreinitDuration (){
    if( interventionMode == forcedEIR ){
        return sim::zero();
    }
    // nYearsWarmupData years for data collection, 50 years stabilization
    return sim::fromYearsI(50) + sim::fromYearsI(nYearsWarmupData);
}
SimTime NonVectorModel::expectedInitDuration (){
    return sim::zero();
}
SimTime NonVectorModel::initIterate (){
    simulationMode = interventionMode;
    if( simulationMode != dynamicEIR ){
        return sim::zero();
    }
    
    // initialKappa is used in calculateEIR
    assert( initialKappa.size() >= sim::stepsPerYear() );
    assert( mod_nn(initialKappa.size(), sim::stepsPerYear()) == 0 );
    for( size_t i=sim::stepsPerYear(); i<initialKappa.size(); ++i ){
        initialKappa[mod_nn(i , sim::stepsPerYear() )] += initialKappa[ i ];
    }
    double factor = static_cast<double>(sim::stepsPerYear()) / static_cast<double>(initialKappa.size());
    initialKappa.resize( sim::stepsPerYear() );
    for(size_t  i = 0; i < initialKappa.size(); ++i) {
        initialKappa[ i ] *= factor;
        // error check:
        if (!(initialKappa[i] > 0.0))     // if not positive
            throw TRACED_EXCEPTION ("initialKappa is invalid", util::Error::InitialKappa);
    }
    
    return sim::zero(); // nothing to do
}

void NonVectorModel::changeEIRIntervention (
        const scnXml::NonVector& nonVectorData)
{
    // Note: requires sim::intervNow() >= sim::zero(), but this can only be
    // called in intervention period anyway.
  simulationMode = transientEIRknown;
  
  if (nSpore != sim::fromDays( nonVectorData.getEipDuration() ))
      throw util::xml_scenario_error ("change-of-EIR intervention cannot change EIP duration");
  
  const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();
  vector<int> nDays( sim::fromDays(daily.size()-1).inSteps() + 1, 0 );
  interventionEIR.assign (nDays.size(), 0.0);
  size_t required_days = static_cast<size_t>(mon::finalSurveyTime().inDays()+1);
  if (daily.size() < required_days) {
    cerr << "Days: " << daily.size() << "\nIntervals: " << nDays.size()
        << "\nRequired: " << required_days << endl;
    throw util::xml_scenario_error ("Insufficient intervention phase EIR values provided");
  }
  //The minimum EIR allowed in the array. The product of the average EIR and a constant.
  double minEIR=min_EIR_mult*averageEIR(nonVectorData);
  for( SimTime mpcday = sim::zero(), endDay = sim::fromDays(daily.size());
      mpcday < endDay; mpcday += sim::oneDay() )
  {
    double EIRdaily = std::max(static_cast<double>(daily[mpcday.inDays()]), minEIR);
    
    // istep is the time period to which the day is assigned.
    size_t istep = mpcday.inSteps();
    nDays[istep]++;
    interventionEIR[istep] += EIRdaily;
  }
  // divide by number of records assigned to each interval (usually one per day)
  for(size_t i = 0; i < interventionEIR.size(); ++i){
    interventionEIR[i] *= sim::oneTS().inDays() / static_cast<double>(nDays[i]);
  }
  
  // I've no idea what this should be, so until someone asks it can be NaN.
  // It was -9.99 and later 0.0. It could of course be recalculated from interventionEIR.
  annualEIR = numeric_limits<double>::quiet_NaN();
}

const map<string,size_t>& NonVectorModel::getSpeciesIndexMap(){
    throw util::xml_scenario_error( "attempt to use a vector-affecting intervention with the non-vector model" );
}

void NonVectorModel::uninfectVectors(){
    if( simulationMode != dynamicEIR )
	cerr <<"Warning: uninfectVectors is not efficacious with forced EIR"<<endl;
    // reset history of human infectivity, which scales dynamic EIR:
    laggedKappa.assign( laggedKappa.size(), 0.0 );
}

void NonVectorModel::deployVectorPopInterv (size_t instance) {
  throw util::xml_scenario_error (viError);
}
void NonVectorModel::deployVectorTrap(size_t instance, double number, SimTime lifespan){
  throw util::xml_scenario_error (viError);
}

void NonVectorModel::update (const Population& population) {
    double currentKappa = TransmissionModel::updateKappa( population );
    
    if( simulationMode == forcedEIR ){
        initialKappa[sim::ts1().moduloSteps(initialKappa.size())] = currentKappa;
    }
}


double NonVectorModel::calculateEIR(Host::Human& human, double ageYears, vector<double>& EIR){
    EIR.resize( 1 );    // no support for per-genotype tracking in this model (possible, but we're lazy)
    // where the full model, with estimates of human mosquito transmission is in use, use this:
    switch (simulationMode) {
        case forcedEIR:
        EIR[0] = initialisationEIR[sim::ts0().moduloYearSteps()];
        break;
        case transientEIRknown:
        // where the EIR for the intervention phase is known, obtain this from
        // the interventionEIR array
        EIR[0] = interventionEIR[sim::intervNow().inSteps()];
        break;
        case dynamicEIR:
        EIR[0] = initialisationEIR[sim::ts0().moduloYearSteps()];
        if (sim::intervNow() >= sim::zero()) {
            // we modulate the initialization based on the human infectiousness time steps ago in the
            // simulation relative to infectiousness at the same time-of-year, pre-intervention.
            // nspore gives the sporozoite development delay.
            size_t t = (sim::ts1()-nSpore).inSteps();
            EIR[0] *=
                laggedKappa[mod_nn(t, laggedKappa.size())] /
                initialKappa[mod_nn(t, sim::stepsPerYear())];
        }
        break;
        default:	// Anything else.. don't continue silently
        throw util::xml_scenario_error ("Invalid simulation mode");
    }
    #ifndef NDEBUG
    if (!(boost::math::isfinite)(EIR[0])) {
        size_t t = (sim::ts1()-nSpore).inSteps();
        ostringstream msg;
        msg << "Error: non-vect eir is: " << EIR[0]
            << "\nlaggedKappa:\t"
            << laggedKappa[mod_nn(t, laggedKappa.size())]
            << "\ninitialKappa:\t"
            << initialKappa[mod_nn(t, sim::stepsPerYear())]
            << endl;
        throw TRACED_EXCEPTION(msg.str(),util::Error::InitialKappa);
    }
    #endif
    EIR[0] *= human.perHostTransmission.relativeAvailabilityHetAge (ageYears);
    return EIR[0];
}


// -----   Private functs ------

double NonVectorModel::averageEIR (const scnXml::NonVector& nonVectorData) {
    // Calculates the arithmetic mean of the whole daily EIR vector read from the .XML file
    double valaverageEIR=0.0;
    size_t i = 0;
    for(const scnXml::NonVector::EIRDailySequence& daily =
        nonVectorData.getEIRDaily(); i < daily.size(); ++i)
    {
        valaverageEIR += (double)daily[i];
    }
    if( i == 0 ) throw util::xml_scenario_error( "no EIRDaily values given" );  // pedantic check
    return valaverageEIR / i;
}


// -----  checkpointing  -----

void NonVectorModel::checkpoint (istream& stream) {
    TransmissionModel::checkpoint (stream);
    nSpore & stream;
    interventionEIR & stream;
    initialKappa & stream;
}
void NonVectorModel::checkpoint (ostream& stream) {
    TransmissionModel::checkpoint (stream);
    nSpore & stream;
    interventionEIR & stream;
    initialKappa & stream;
}

} }
