/*

 This file is part of OpenMalaria.
 
 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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
#include "Transmission/TransmissionModel.h"
#include "Transmission/NonVector.h"
#include "Transmission/Vector/VectorTransmission.h"
#include "Transmission/PerHostTransmission.h"

#include "inputData.h"
#include "Simulation.h"
#include "util/BoincWrapper.h"

#include <cmath> 
#include <cfloat>
#include <gsl/gsl_vector.h>

namespace OM { namespace Transmission {
    
TransmissionModel* TransmissionModel::createTransmissionModel () {
  // EntoData contains either a list of at least one anopheles or a list of at
  // least one EIRDaily.
  const scnXml::EntoData::VectorOptional& vectorData = InputData().getEntoData().getVector();
  if (vectorData.present())
    return new VectorTransmission(vectorData.get());
  else {
      const scnXml::EntoData::NonVectorOptional& nonVectorData = InputData().getEntoData().getNonVector();
    if (!nonVectorData.present())	// should be a validation error, but anyway...
      throw util::xml_scenario_error ("Neither vector nor non-vector data present in the XML!");
    return new NonVectorTransmission(nonVectorData.get());
  }
}

TransmissionModel::TransmissionModel() :
    ageCorrectionFactor(numeric_limits<double>::signaling_NaN()),
    simulationMode(equilibriumMode),
    _sumAnnualKappa(0.0), annualEIR(0.0),
    timeStepNumEntoInnocs (0)
#ifdef OMV_CSV_REPORTING
  , csvReporting ("vector.txt", ios::app)
#endif
{
  kappa.resize (Global::intervalsPerYear, 0.0);
  initialisationEIR.resize (Global::intervalsPerYear);
  innoculationsPerAgeGroup.resize (SurveyAgeGroup::getNumGroups(), 0.0);
  innoculationsPerDayOfYear.resize (Global::intervalsPerYear, 0.0);
  timeStepEntoInnocs.resize (SurveyAgeGroup::getNumGroups(), 0.0);
  
  // noOfAgeGroupsSharedMem must be at least as large as both of these to avoid
  // memory corruption or extra tests when setting/copying values
  noOfAgeGroupsSharedMem = std::max(SurveyAgeGroup::getNumGroups(), util::SharedGraphics::KappaArraySize);
}

TransmissionModel::~TransmissionModel () {
#ifdef OMV_CSV_REPORTING
  csvReporting.close();
#endif
}


void TransmissionModel::updateAgeCorrectionFactor (const std::list<Host::Human>& population, int populationSize) {
    // Calculate relative availability correction, so calls from vectorUpdate,
    // etc., will have a mean of 1.0.
    double meanRelativeAvailability = 0.0;
    for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h)
	meanRelativeAvailability += Transmission::PerHostTransmission::relativeAvailabilityAge (h->getAgeInYears());
    ageCorrectionFactor = populationSize / meanRelativeAvailability;
}

void TransmissionModel::updateKappa (const std::list<Host::Human>& population, int simulationTime) {
  // We calculate kappa for output and non-vector model, and kappaByAge for
  // the shared graphics.
  
  double sumWt_kappa= 0.0;
  double sumWeight  = 0.0;
  kappaByAge.assign (noOfAgeGroupsSharedMem, 0.0);
  nByAge.assign (noOfAgeGroupsSharedMem, 0);
  
  for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
    double ageYears = h->getAgeInYears();
    double t = h->perHostTransmission.relativeAvailabilityHetAge(ageYears) * ageCorrectionFactor;
    sumWeight += t;
    t *= h->probTransmissionToMosquito();
    sumWt_kappa += t;
    
    // kappaByAge and nByAge are used in the screensaver only
    SurveyAgeGroup ag = h->ageGroup();
    kappaByAge[ag.i()] += t;
    ++nByAge[ag.i()];
  }
  
#ifndef NDEBUG
  if (sumWeight < DBL_MIN * 4.0)	// if approx. eq. 0 or negative
    throw runtime_error ("sumWeight is invalid");
#endif
  
  size_t tmod = (simulationTime-1) % Global::intervalsPerYear;
  kappa[tmod] = sumWt_kappa / sumWeight;
  
  //Calculate time-weighted average of kappa
  _sumAnnualKappa += kappa[tmod] * initialisationEIR[tmod];
  if (tmod == Global::intervalsPerYear - 1) {
    if (annualEIR == 0.0) {
      _annualAverageKappa=0;
      cerr << "aE.eq.0" << endl;
    } else {
      _annualAverageKappa = _sumAnnualKappa / annualEIR;
    }
    _sumAnnualKappa = 0.0;
  }
  
  double timeStepTotal = 0.0;
  for (size_t group = 0; group < timeStepEntoInnocs.size(); ++group) {
    timeStepTotal += timeStepEntoInnocs[group];
    innoculationsPerAgeGroup[group] += timeStepEntoInnocs[group];
    // Reset to zero:
    timeStepEntoInnocs[group] = 0.0;
  }
  innoculationsPerDayOfYear[tmod] = timeStepTotal / timeStepNumEntoInnocs;
  timeStepNumEntoInnocs = 0;
  
  // Shared graphics: report infectiousness
  if (Global::simulationTime % 6 ==  0) {
    for (size_t i = 0; i < noOfAgeGroupsSharedMem; i++)
      kappaByAge[i] /= nByAge[i];
    util::SharedGraphics::copyKappa(&kappaByAge[0]);
  }
  
#ifdef OMV_CSV_REPORTING
  csvReporting << initialisationEIR[simulationTime%Global::intervalsPerYear]
	       << '\t' << innoculationsPerDayOfYear[tmod]
	       << '\t' << kappa[tmod]
	       << '\t' << ageCorrectionFactor
	       << '\t'<< endl;
#endif
}

double TransmissionModel::getEIR (int simulationTime, PerHostTransmission& host, double ageInYears, SurveyAgeGroup ageGroup) {
  /* For the NonVector model, the EIR should just be multiplied by the
   * availability. For the Vector model, the availability is also required
   * for internal calculations, but again the EIR should be multiplied by the
   * availability. */
  double EIR = calculateEIR (simulationTime, host, ageInYears);
  
  timeStepEntoInnocs[ageGroup.i()] += EIR;
  timeStepNumEntoInnocs ++;
  return EIR;
}

void TransmissionModel::summarize (Survey& survey) {
  survey.setNumTransmittingHosts(kappa[(Global::simulationTime-1) % Global::intervalsPerYear]);
  survey.setAnnualAverageKappa(_annualAverageKappa);
  
  survey.setInnoculationsPerAgeGroup (innoculationsPerAgeGroup);	// Array contents must be copied.
  innoculationsPerAgeGroup.assign (innoculationsPerAgeGroup.size(), 0.0);
  
  survey.set_Vector_EIR_Input (initialisationEIR[Global::simulationTime % Global::intervalsPerYear]);
  survey.set_Vector_EIR_Simulated (innoculationsPerDayOfYear[(Global::simulationTime-1) % Global::intervalsPerYear]);
}

void TransmissionModel::intervLarviciding (const scnXml::Larviciding&) {
  throw util::xml_scenario_error ("larviciding when not using a Vector model");
}


// -----  checkpointing  -----

void TransmissionModel::checkpoint (istream& stream) {
    ageCorrectionFactor & stream;
    simulationMode & stream;
    initialisationEIR & stream;
    kappa & stream;
    _annualAverageKappa & stream;
    _sumAnnualKappa & stream;
    annualEIR & stream;
    innoculationsPerDayOfYear & stream;
    innoculationsPerAgeGroup & stream;
    timeStepEntoInnocs & stream;
    timeStepNumEntoInnocs & stream;
    noOfAgeGroupsSharedMem & stream;
    kappaByAge & stream;
    nByAge & stream;
}
void TransmissionModel::checkpoint (ostream& stream) {
    ageCorrectionFactor & stream;
    simulationMode & stream;
    initialisationEIR & stream;
    kappa & stream;
    _annualAverageKappa & stream;
    _sumAnnualKappa & stream;
    annualEIR & stream;
    innoculationsPerDayOfYear & stream;
    innoculationsPerAgeGroup & stream;
    timeStepEntoInnocs & stream;
    timeStepNumEntoInnocs & stream;
    noOfAgeGroupsSharedMem & stream;
    kappaByAge & stream;
    nByAge & stream;
}

} }