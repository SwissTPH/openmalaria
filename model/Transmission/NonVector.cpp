/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#include "Transmission/NonVector.h"
#include "Transmission/PerHostTransmission.h"
#include "inputData.h"
#include "util/random.h"
#include "util/vectors.h"
#include "util/StreamValidator.h"
#include "Monitoring/Surveys.h"	// sim-end timestep
#include <limits>
#include <cmath>

namespace OM { namespace Transmission {
    namespace vectors = util::vectors;

//static (class) variables
const double NonVectorTransmission::totalInfectionrateVariance= 1.0;
const double NonVectorTransmission::min_EIR_mult= 0.01;

const int nYearsWarmupData = 5;

NonVectorTransmission::NonVectorTransmission(const scnXml::NonVector& nonVectorData) :
  nspore( TimeStep::fromDays( nonVectorData.getEipDuration() ) )
{
    laggedKappa.resize( nspore.asInt()+1, 0.0 );
    
    vector<int> nDays (TimeStep::stepsPerYear, 0);
    //The minimum EIR allowed in the array. The product of the average EIR and a constant.
    double minEIR=min_EIR_mult*averageEIR(nonVectorData);
    
    const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();
    if( daily.size() < static_cast<size_t>(TimeStep::intervalsPerYear.inDays()) )
        throw util::xml_scenario_error( "insufficient EIRDaily data for a year" );
    
    for (size_t mpcday = 0; mpcday < daily.size(); ++mpcday) {
        double EIRdaily = std::max((double)daily[mpcday], minEIR);
        
        // istep is the time period to which the day is assigned.  The result of the
        // division is automatically rounded down to the next integer.
        size_t i1 = (mpcday / TimeStep::interval) % TimeStep::stepsPerYear;
        //EIR() is the sum of the EIRs assigned to the 73 different recurring time points
        nDays[i1]++;
        initialisationEIR[i1] += EIRdaily;
    }
    
    // Calculate total annual EIR
    // divide by number of records assigned to each interval (usually one per day)
    for (TimeStep j(0);j<TimeStep::intervalsPerYear; ++j) {
        initialisationEIR[j.asInt()] *= TimeStep::interval / (double)nDays[j.asInt()];
        annualEIR += initialisationEIR[j.asInt()];
    }
    
    initialKappa.assign( TimeStep::fromYears(nYearsWarmupData).asInt(), 0.0 );
}

NonVectorTransmission::~NonVectorTransmission () {}

void NonVectorTransmission::scaleEIR (double factor){
    vectors::scale( initialisationEIR, factor );
    annualEIR = vectors::sum( initialisationEIR );
}
void NonVectorTransmission::scaleXML_EIR (scnXml::EntoData& ed, double factor) const{
    assert( ed.getNonVector().present() );
    scnXml::NonVector::EIRDailySequence& daily = ed.getNonVector().get().getEIRDaily();
    
    for( scnXml::NonVector::EIRDailyIterator it = daily.begin();
	it != daily.end(); ++it ){
	double old = *it;
	*it = old * factor;
    }
}


TimeStep NonVectorTransmission::minPreinitDuration (){
    if( InputData().getEntomology().getMode() == equilibriumMode ){
        return TimeStep(0);
    }
    // nYearsWarmupData years for data collection, 50 years stabilization
    return TimeStep::fromYears(50) + TimeStep::fromYears(nYearsWarmupData);
}
TimeStep NonVectorTransmission::expectedInitDuration (){
    return TimeStep(0);
}
TimeStep NonVectorTransmission::initIterate (){
    simulationMode = interventionMode;
    if( simulationMode != dynamicEIR ){
        return TimeStep(0);
    }
    
    // initialKappa is used in calculateEIR
    size_t yearLen = TimeStep::fromYears(1).asInt();
    assert( initialKappa.size() >= yearLen );
    assert( initialKappa.size() % yearLen == 0 );
    for( size_t i=yearLen; i<initialKappa.size(); ++i ){
        initialKappa[ i % yearLen ] += initialKappa[ i ];
    }
    double factor = static_cast<double>(yearLen) / static_cast<double>(initialKappa.size());
    initialKappa.resize( yearLen );
    for (size_t  i = 0; i < initialKappa.size(); ++i) {
        initialKappa[ i ] *= factor;
        // error check:
        if (!(initialKappa[i] > 0.0))     // if not positive
            throw util::traced_exception ("initialKappa is invalid");
    }
    
    return TimeStep(0); // nothing to do
}

void NonVectorTransmission::setTransientEIR (const scnXml::NonVector& nonVectorData) {
    // Note: requires TimeStep::interventionPeriod >= 0, but this can only be called in intervention period anyway.
  simulationMode = transientEIRknown;
  
  if (nspore != TimeStep::fromDays( nonVectorData.getEipDuration() ))
      throw util::xml_scenario_error ("change-of-EIR intervention cannot change EIP duration");
  
  const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();
  vector<int> nDays ((daily.size()-1)/TimeStep::interval + 1, 0);
  interventionEIR.assign (nDays.size(), 0.0);
  if (daily.size() < static_cast<size_t>(Monitoring::Surveys.getFinalTimestep().inDays()+1)) {
    cerr << "Days: " << daily.size() << "\nIntervals: " << nDays.size() << "\nRequired: " << Monitoring::Surveys.getFinalTimestep().inDays()+1 << endl;
    throw util::xml_scenario_error ("Insufficient intervention phase EIR values provided");
  }
  //The minimum EIR allowed in the array. The product of the average EIR and a constant.
  double minEIR=min_EIR_mult*averageEIR(nonVectorData);
  for (size_t mpcday = 0; mpcday < daily.size(); ++mpcday) {
    double EIRdaily = std::max((double)daily[mpcday], minEIR);
    
    // istep is the time period to which the day is assigned.  The result of the
    // division is automatically rounded down to the next integer.
    size_t istep = mpcday / TimeStep::interval;
    nDays[istep]++;
    interventionEIR[istep] += EIRdaily;
  }
  // divide by number of records assigned to each interval (usually one per day)
  for (size_t i = 0; i < interventionEIR.size(); ++i)
    interventionEIR[i] *= TimeStep::interval / nDays[i];
  
  // I've no idea what this should be, so until someone asks it can be NaN.
  // It was -9.99 and later 0.0. It could of course be recalculated from interventionEIR.
  annualEIR = numeric_limits<double>::quiet_NaN();
}

void NonVectorTransmission::changeEIRIntervention (const scnXml::NonVector& ed) {
    setTransientEIR (ed);
}

void NonVectorTransmission::uninfectVectors(){
    if( simulationMode != dynamicEIR )
	cerr <<"Warning: uninfectVectors is not efficacious with forced EIR"<<endl;
    // reset history of human infectivity, which scales dynamic EIR:
    laggedKappa.assign( laggedKappa.size(), 0.0 );
}

void NonVectorTransmission::update (const std::list<Host::Human>& population, int populationSize) {
    double currentKappa = TransmissionModel::updateKappa( population );
    
    if( simulationMode == equilibriumMode ){
        initialKappa[ TimeStep::simulation % initialKappa.size() ] = currentKappa;
    }
}


double NonVectorTransmission::calculateEIR(PerHostTransmission& perHost, double ageYears){
  // where the full model, with estimates of human mosquito transmission is in use, use this:
  double eir;
  switch (simulationMode) {
    case equilibriumMode:
      //FIXME: figure out if it is correct that 1 is subtracted here and below but not elsewhere
      eir = initialisationEIR[(TimeStep::simulation-TimeStep(1)) % TimeStep::stepsPerYear];
      break;
    case transientEIRknown:
      // where the EIR for the intervention phase is known, obtain this from
      // the interventionEIR array (why -1? See interventionEIR declaration)
      eir = interventionEIR[TimeStep::interventionPeriod.asInt() - 1];
      break;
    case dynamicEIR:
      eir = initialisationEIR[(TimeStep::simulation-TimeStep(1)) % TimeStep::stepsPerYear];
      if (TimeStep::interventionPeriod >= TimeStep(0)) {
	  // we modulate the initialization based on the human infectiousness  timesteps ago in the
	  // simulation relative to infectiousness at the same time-of-year, pre-intervention.
	  // nspore gives the sporozoite development delay.
	eir *=
            laggedKappa[(TimeStep::simulation-nspore) % laggedKappa.size()] /
            initialKappa[(TimeStep::simulation-nspore) % TimeStep::stepsPerYear];
      }
      break;
    default:	// Anything else.. don't continue silently
      throw util::xml_scenario_error ("Invalid simulation mode");
  }
#ifndef NDEBUG
  if (!finite(eir)) {
    ostringstream msg;
    msg << "Error: non-vect eir is: " << eir
	<< "\nlaggedKappa:\t" << laggedKappa[(TimeStep::simulation-nspore) % laggedKappa.size()]
	<< "\ninitialKappa:\t" << initialKappa[(TimeStep::simulation-nspore) % TimeStep::stepsPerYear] << endl;
    throw util::traced_exception(msg.str(),util::Error::InitialKappa);
  }
#endif
  return eir * perHost.relativeAvailabilityHetAge (ageYears);
}


// -----   Private functs ------

double NonVectorTransmission::averageEIR (const scnXml::NonVector& nonVectorData) {
  // Calculates the arithmetic mean of the whole daily EIR vector read from the .XML file
  double valaverageEIR=0.0;
  size_t i = 0;
  for (const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();
       i < daily.size(); ++i) {
    valaverageEIR += (double)daily[i];
  }
  return valaverageEIR / i;
}


// -----  checkpointing  -----

void NonVectorTransmission::checkpoint (istream& stream) {
    TransmissionModel::checkpoint (stream);
    nspore & stream;
    interventionEIR & stream;
    initialKappa & stream;
}
void NonVectorTransmission::checkpoint (ostream& stream) {
    TransmissionModel::checkpoint (stream);
    nspore & stream;
    interventionEIR & stream;
    initialKappa & stream;
}

} }
