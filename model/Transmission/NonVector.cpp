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
#include "Transmission/PerHost.h"
#include "inputData.h"
#include "simulation.h"
#include "GSLWrapper.h"
#include <cfloat>

//static (class) variables
const double NonVectorTransmission::totalInfectionrateVariance= 1.0;
const double NonVectorTransmission::min_EIR_mult= 0.01; 

NonVectorTransmission::NonVectorTransmission(const scnXml::NonVector& nonVectorData)
{
  nspore = nonVectorData.getEipDuration() / Global::interval;
  
  initialKappa.resize (Global::intervalsPerYear);
  
  vector<int> nDays (Global::intervalsPerYear, 0);
  initialisationEIR.assign (Global::intervalsPerYear, 0.0);
  //The minimum EIR allowed in the array. The product of the average EIR and a constant.
  double minEIR=min_EIR_mult*averageEIR(nonVectorData);
  const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();
  for (size_t mpcday = 0; mpcday < daily.size(); ++mpcday) {
    double EIRdaily = std::max((double)daily[mpcday], minEIR);
    
    // istep is the time period to which the day is assigned.  The result of the
    // division is automatically rounded down to the next integer.
    size_t i1 = (mpcday / Global::interval) % Global::intervalsPerYear;
    //EIR() is the sum of the EIRs assigned to the 73 different recurring time points
    nDays[i1]++;
    initialisationEIR[i1] += EIRdaily;
  }
  
  // Calculate total annual EIR
  for (size_t j=0;j<Global::intervalsPerYear; j++) {
    // NOTE: this should be unnecessary, if daily.size() == daysInYear
    initialisationEIR[j] *= Global::interval / (double)nDays[j];
    annualEIR += initialisationEIR[j];
  }
}

NonVectorTransmission::~NonVectorTransmission () {}


//! initialise the main simulation 
void NonVectorTransmission::initMainSimulation (const std::list<Human>&, int){
  // initialKappa is used in calculateEIR
  copyToInitialKappa();
}


void NonVectorTransmission::setTransientEIR (const scnXml::NonVector& nonVectorData) {
  Global::simulationMode = transientEIRknown;
  const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();
  vector<int> nDays ((daily.size()-1)/Global::interval + 1, 0);
  interventionEIR.assign (nDays.size(), 0.0);
  if (nDays.size() <= Simulation::simulationDuration) {
    cerr << "Days: " << daily.size() << "\nIntervals: " << nDays.size() << "\nRequired: " << Simulation::simulationDuration+1 << endl;
    throw xml_scenario_error ("Insufficient intervention phase EIR values provided");
  }
  //The minimum EIR allowed in the array. The product of the average EIR and a constant.
  double minEIR=min_EIR_mult*averageEIR(nonVectorData);
  for (size_t mpcday = 0; mpcday < daily.size(); ++mpcday) {
    double EIRdaily = std::max((double)daily[mpcday], minEIR);
    
    // istep is the time period to which the day is assigned.  The result of the
    // division is automatically rounded down to the next integer.
    size_t istep = mpcday / Global::interval;
    nDays[istep]++;
    interventionEIR[istep] += EIRdaily;
  }
  // NOTE: this should be unnecessary, if daily.size() == daysInYear
  for (size_t i = 0; i < interventionEIR.size(); ++i)
    interventionEIR[i] *= Global::interval / nDays[i];
  annualEIR=-9.99;
}

void NonVectorTransmission::copyToInitialKappa () {
  initialKappa = kappa;
  //NOTE: error check
  for (size_t  i = 0; i < initialKappa.size(); ++i) {
    if (initialKappa[i] < DBL_MIN * 4.0)	// if approx. eq. 0 or negative
      throw range_error ("initialKappa is invalid");
  }
}

double NonVectorTransmission::calculateEIR(int simulationTime, PerHostTransmission& perHost, double ageInYears){
  // where the full model, with estimates of human mosquito transmission is in use, use this:
  double eir;
  switch (Global::simulationMode) {
    case transientEIRknown:
      // where the EIR for the intervention phase is known, obtain this from
      // the interventionEIR array
      eir = interventionEIR[Simulation::timeStep];
      break;
    case dynamicEIR:
      eir = initialisationEIR[(simulationTime-1) % Global::intervalsPerYear];
      if (Simulation::timeStep != 0) {
	eir *=
            kappa[(simulationTime-nspore-1) % Global::intervalsPerYear] /
            initialKappa[(simulationTime-nspore-1) % Global::intervalsPerYear];
      }
      break;
    default:	// Anything else.. don't continue silently
      throw xml_scenario_error ("Invalid simulation mode");
  }
#ifndef NDEBUG
  if (!finite(eir)) {
    ostringstream msg;
    msg << "Error: non-vect eir is: " << eir
	<< "\nkappa:\t" << kappa[(simulationTime-nspore-1) % Global::intervalsPerYear]
	<< "\ninitialKappa:\t" << initialKappa[(simulationTime-nspore-1) % Global::intervalsPerYear] << endl;
    throw overflow_error(msg.str());
  }
#endif
  return eir * getRelativeAvailability(ageInYears) * perHost.entoAvailability();
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
