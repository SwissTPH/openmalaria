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
#include "TransmissionModel/NonVector.h"
#include "TransmissionModel/PerHost.h"
#include "inputData.h"
#include "simulation.h"
#include "GSLWrapper.h"

//static (class) variables
const double NonVectorTransmission::totalInfectionrateVariance= 1.0;
const double NonVectorTransmission::min_EIR_mult= 0.01; 

NonVectorTransmission::NonVectorTransmission(const scnXml::NonVector& nonVectorData)
{
  nspore = nonVectorData.getEipDuration() / Global::interval;
  maxIntervals=maxDurIntPhaseEIR / Global::interval;
  
  initialKappa.resize (Global::intervalsPerYear);
  inputEIR(nonVectorData);
}

NonVectorTransmission::~NonVectorTransmission () {}


//! initialise the main simulation 
void NonVectorTransmission::initMainSimulation (int populationSize){
  // initialKappa is used in calculateEIR
  copyToInitialKappa();
}


void NonVectorTransmission::inputEIR (const scnXml::NonVector& nonVectorData) {
  vector<int> nDays;
  //initialise all the EIR arrays to 0
  if (Global::simulationMode != transientEIRknown) {
    initialisationEIR.assign (Global::intervalsPerYear, 0.0);
    nDays.assign (Global::intervalsPerYear, 0);
  } else {
    intEIR.assign (maxIntervals, 0.0);
    nDays.assign (maxIntervals, 0);
  }
  //The minimum EIR allowed in the array. The product of the average EIR and a constant.
  double minEIR=min_EIR_mult*averageEIR(nonVectorData);
  const scnXml::NonVector::EIRDailySequence& daily = nonVectorData.getEIRDaily();
  for (size_t mpcday = 0; mpcday < daily.size(); ++mpcday) {
    // FIXME: convertions to make this run the same as when encapsulated in updateEIR
    int day = mpcday;
    double EIRdaily = std::max((double)daily[mpcday], minEIR);
    
    // istep is the time period to which the day is assigned.  The result of the
    // division is automatically rounded down to the next integer.
    // Except: why is 1 subtracted from day? -1/5 is usually be rounded to 0, but
    // some compilers may round it to -1.
    int istep = (day-1) / Global::interval;
    if (Global::simulationMode !=  transientEIRknown) {
      int i1 = Global::modIntervalsPerYear(1+istep) - 1;
      size_t i2 = (1+istep) % Global::intervalsPerYear;
      nDays[i1]++;
      //EIR() is the arithmetic mean of the EIRs assigned to the 73 different recurring time points
      initialisationEIR[i1] = ((initialisationEIR[i1] * (nDays[i1]-1)) + EIRdaily) / nDays[i1];
    } else {
      nDays[istep]++;
      intEIR[istep]= ((intEIR[istep] * (nDays[istep]-1)) + EIRdaily) / nDays[istep];
    }
  }
  
  // Calculate total annual EIR
  if (Global::simulationMode != transientEIRknown) {
    annualEIR=0.0;
    for (size_t j=0;j<Global::intervalsPerYear; j++) {
      annualEIR += Global::interval*initialisationEIR[j];
    }
  } else {
    annualEIR=-9.99;
  }
}

void NonVectorTransmission::copyToInitialKappa () {
  initialKappa = kappa;
}

double NonVectorTransmission::calculateEIR(int simulationTime, PerHostTransmission&){
  // where the full model, with estimates of human mosquito transmission is in use, use this:
  switch (Global::simulationMode) {
    case transientEIRknown:
      // where the EIR for the intervention phase is known, obtain this from
      // the intEIR array
      return intEIR[Simulation::timeStep - 1];
      break;
    case dynamicEIR:
      if (Simulation::timeStep == 1) {
	return initialisationEIR[(simulationTime-1) % Global::intervalsPerYear];
      } else {
	return initialisationEIR[(simulationTime-1) % Global::intervalsPerYear] *
            kappa[(simulationTime-nspore-1) % Global::intervalsPerYear] /
            initialKappa[(simulationTime-nspore-1) % Global::intervalsPerYear];
      }
      break;
    default:	// Anything else.. don't continue silently
      throw xml_scenario_error ("Invalid simulation mode");
  }
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
