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
#include "intervention.h"

//static (class) variables
const double NonVectorTransmission::susceptibility= 0.702;
const double NonVectorTransmission::totalInfectionrateVariance= 1.0;
const double NonVectorTransmission::min_EIR_mult= 0.01; 

NonVectorTransmission::NonVectorTransmission(const scnXml::NonVector& nonVectorData)
{
  nspore = nonVectorData.getEipDuration() / Global::interval;
  gamma_p=getParameter(Params::GAMMA_P);
  Sinf=1-exp(-getParameter(Params::NEG_LOG_ONE_MINUS_SINF));
  Simm=getParameter(Params::SIMM);
  Estar=getParameter(Params::E_STAR);
  Xstar_p=getParameter(Params::X_STAR_P);
  nDays = (int *) malloc(((Global::intervalsPerYear))*sizeof(int));
  maxIntervals=maxDurIntPhaseEIR / Global::interval;
  ino = (int *) malloc(((maxIntervals))*sizeof(int));
  intEIR = (double *) malloc(((maxIntervals))*sizeof(double));
  
  //! constant defining the constraint for the Gamma shape parameters
  /// Used for the case where availability is assumed gamma distributed
  double r_square_Gamma;
  /*
  //! Expected number of inoculations
  /// Product of measured EIR, susceptibility and length of time Global::interval
  double gsi = 1.0;
  
  r_square_Gamma=(totalInfectionrateVariance**2-gsi*BaselineAvailabilityMean)/(gsi*BaselineAvailabilityMean)**2
  r_square_Gamma must be greater than zero, so r_square_LogNormal is also. 
  */
  r_square_Gamma=0.649;
  //such that r_square_LogNormal =0.5
  
  //! constant defining the constraint for the log Normal variance
  /// Used for the case where availability is assumed log Normally distributed
  double r_square_LogNormal = log(1.0+r_square_Gamma);
  
  //TODO: Sanity check for sqrt and division by zero
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    InfectionrateShapeParam = (PerHostTransmission::BaselineAvailabilityShapeParam+1.0) / (r_square_Gamma*PerHostTransmission::BaselineAvailabilityShapeParam - 1.0);
    InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  }
  else if (Global::modelVersion &
           (LOGNORMAL_MASS_ACTION | LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM)) {
    InfectionrateShapeParam = sqrt(r_square_LogNormal - 1.86*pow(PerHostTransmission::BaselineAvailabilityShapeParam, 2));
    InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  }
  
  inputEIR(nonVectorData);
}

NonVectorTransmission::~NonVectorTransmission () {
  free(nDays);
  free(ino);
  free(intEIR);
}


//! initialise the main simulation 
void NonVectorTransmission::initMainSimulation (int populationSize){
  // initialKappa is used in calculateEIR
  copyToInitialKappa();
}


void NonVectorTransmission::inputEIR (const scnXml::NonVector& nonVectorData) {
  //initialise all the EIR arrays to 0
  if (Global::simulationMode != transientEIRknown) {
    for (size_t j=0;j<Global::intervalsPerYear; j++) {
      initialisationEIR[j]=0.0;
      nDays[j]=0;
    }
  } else {
    for (int j=0;j<maxIntervals; j++) {
      intEIR[j]=0.0;
      ino[j]=0;
    }
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
      size_t i2 = (1+istep) % Global::intervalsPerYear;
      nDays[i2]++;
      //EIR() is the arithmetic mean of the EIRs assigned to the 73 different recurring time points
      initialisationEIR[i2] = ((initialisationEIR[i2] * (nDays[i2]-1)) + EIRdaily) / nDays[i2];
    } else {
      ino[istep]++;
      intEIR[istep]= ((intEIR[istep] * (ino[istep]-1)) + EIRdaily) / ino[istep];
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
	return initialisationEIR[simulationTime % Global::intervalsPerYear];
      } else {
	return initialisationEIR[simulationTime % Global::intervalsPerYear] *
            kappa[(simulationTime-nspore) % Global::intervalsPerYear] /
            initialKappa[(simulationTime-nspore) % Global::intervalsPerYear];
      }
      break;
    default:	// Anything else.. don't continue silently
      throw xml_scenario_error ("Invalid simulation mode");
  }
}


double NonVectorTransmission::getExpectedNumberOfInfections (Human& human, double age_adj_EIR) {
  PerHostTransmission& host = human._perHostTransmission;
  double baseAvailToMos = host._BaselineAvailabilityToMosquitoes;
  //The age-adjusted EIR, possibly adjusted for bed nets.
  double expectedNumInfections;
  
  double ExpectedInfectionRate = age_adj_EIR * baseAvailToMos * susceptibility * Global::interval;
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    expectedNumInfections = (W_GAMMA((InfectionrateShapeParam), (ExpectedInfectionRate/InfectionrateShapeParam)));
  } else if (Global::modelVersion & LOGNORMAL_MASS_ACTION) {
    expectedNumInfections = sampleFromLogNormal(W_UNIFORM(),
        log(ExpectedInfectionRate) - 0.5*pow(InfectionrateShapeParam, 2),
        InfectionrateShapeParam);
  } else {
    //The default model is that in Smith et al, AJTMH 2006 75 Suppl 2
    double survivalOfInoculum=(1.0+pow((host._cumulativeEIRa/Xstar_p), gamma_p));
    survivalOfInoculum = Simm+(1.0-Simm)/survivalOfInoculum;
    survivalOfInoculum = survivalOfInoculum*(Sinf+(1-Sinf)/(1 + age_adj_EIR/Estar));
    
    if(Global::modelVersion & LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM) {
      expectedNumInfections = survivalOfInoculum *
          sampleFromLogNormal(W_UNIFORM(),
                              log(ExpectedInfectionRate) - 0.5*pow(InfectionrateShapeParam, 2),
                              InfectionrateShapeParam);
    } else {
      expectedNumInfections = survivalOfInoculum *
          age_adj_EIR * Global::interval * baseAvailToMos;
    }
  }
  
  //Introduce the effect of vaccination. Note that this does not affect cumEIR.
  if (Vaccine::PEV.active) {
    expectedNumInfections *= (1 - human.getPEVEfficacy());
  }
  return expectedNumInfections;
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
