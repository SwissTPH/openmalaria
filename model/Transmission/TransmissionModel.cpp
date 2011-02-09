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
#include "Monitoring/Continuous.h"
#include "util/BoincWrapper.h"
#include "util/StreamValidator.h"
#include "util/CommandLine.h"
#include "util/vectors.h"

#include <cmath> 
#include <cfloat>
#include <gsl/gsl_vector.h>

namespace OM { namespace Transmission {
namespace vectors = util::vectors;

TransmissionModel* TransmissionModel::createTransmissionModel (int populationSize) {
  // EntoData contains either a list of at least one anopheles or a list of at
  // least one EIRDaily.
  const scnXml::EntoData& entoData = InputData().getEntoData();
  const scnXml::EntoData::VectorOptional& vectorData = entoData.getVector();
  
  TransmissionModel *model;
  if (vectorData.present())
    model = new VectorTransmission(vectorData.get(), populationSize);
  else {
      const scnXml::EntoData::NonVectorOptional& nonVectorData = entoData.getNonVector();
    if (!nonVectorData.present())       // should be a validation error, but anyway...
      throw util::xml_scenario_error ("Neither vector nor non-vector data present in the XML!");
    model = new NonVectorTransmission(nonVectorData.get());
  }
  
  if( entoData.getAnnualEIR().present() ){
      model->scaleEIR( entoData.getAnnualEIR().get() / model->annualEIR );
      assert( vectors::approxEqual( model->annualEIR, entoData.getAnnualEIR().get() ) );
  }
  
  if( util::CommandLine::option( util::CommandLine::PRINT_ANNUAL_EIR ) ){
      //Note: after internal scaling (which doesn't imply exit)
      //but before external scaling.
      cout << "Total annual EIR: "<<model->annualEIR<<endl;
  }
  if( util::CommandLine::option( util::CommandLine::SET_ANNUAL_EIR ) ){
      model->scaleXML_EIR(
        InputData.getMutableScenario().getEntoData(),
        util::CommandLine::getNewEIR() / model->annualEIR
      );
      InputData.documentChanged = true;
  }
  
  return model;
}


// The times here should be for the last updated index of arrays:
void TransmissionModel::ctsCbInputEIR (ostream& stream){
    stream<<'\t'<<initialisationEIR[TimeStep::simulation % TimeStep::stepsPerYear];
}
void TransmissionModel::ctsCbSimulatedEIR (ostream& stream){
    stream<<'\t'<<innoculationsPerDayOfYear[TimeStep::simulation % TimeStep::stepsPerYear];
}
void TransmissionModel::ctsCbKappa (ostream& stream){
    stream<<'\t'<<kappa[TimeStep::simulation % TimeStep::stepsPerYear];
}
void TransmissionModel::ctsCbHumanAvail (ostream& stream){
    stream<<'\t'<<1.0/ageCorrectionFactor;
}
void TransmissionModel::ctsCbNumTransmittingHumans (ostream& stream){
    stream<<'\t'<<numTransmittingHumans;
}

TransmissionModel::TransmissionModel() :
    ageCorrectionFactor(numeric_limits<double>::signaling_NaN()),
    simulationMode(equilibriumMode),
    interventionMode(InputData().getEntoData().getMode()),
    _sumAnnualKappa(0.0),
    BSSInitialisationEIR(0.0), BSSInnoculationsPerDayOfYear(0.0), BSSTimesteps(0),
    annualEIR(0.0),
    timeStepNumEntoInnocs (0)
{
    if (interventionMode != equilibriumMode && interventionMode != dynamicEIR){
        // Note: previously 3 was allowed -- but mode is set to 3 anyway when
        // "intervention" EIR data is loaded, so 2 or 4 should be used here.
        throw util::xml_scenario_error("mode attribute has invalid value (expected: 2 or 4)");
    }
    
  kappa.assign (TimeStep::stepsPerYear, 0.0);
  initialisationEIR.assign (TimeStep::stepsPerYear, 0.0);
  innoculationsPerAgeGroup.assign (Monitoring::AgeGroup::getNumGroups(), 0.0);
  innoculationsPerDayOfYear.assign (TimeStep::stepsPerYear, 0.0);
  timeStepEntoInnocs.assign (Monitoring::AgeGroup::getNumGroups(), 0.0);
  
  // noOfAgeGroupsSharedMem must be at least as large as both of these to avoid
  // memory corruption or extra tests when setting/copying values
  noOfAgeGroupsSharedMem = std::max(Monitoring::AgeGroup::getNumGroups(), util::SharedGraphics::KappaArraySize);
  
  using Monitoring::Continuous;
  Continuous::registerCallback( "input EIR", "\tinput EIR", MakeDelegate( this, &TransmissionModel::ctsCbInputEIR ) );
  Continuous::registerCallback( "simulated EIR", "\tsimulated EIR", MakeDelegate( this, &TransmissionModel::ctsCbSimulatedEIR ) );
  Continuous::registerCallback( "human infectiousness", "\thuman infectiousness", MakeDelegate( this, &TransmissionModel::ctsCbKappa ) );
  Continuous::registerCallback( "human availability", "\tmean human availability", MakeDelegate( this, &TransmissionModel::ctsCbHumanAvail ) );
  Continuous::registerCallback( "num transmitting humans", "\tnum transmitting humans", MakeDelegate( this, &TransmissionModel::ctsCbNumTransmittingHumans ) );
}

TransmissionModel::~TransmissionModel () {
}


void TransmissionModel::updateAgeCorrectionFactor (std::list<Host::Human>& population, int populationSize) {
    //NOTE: ageCorrectionFactor is now only used within the vector
    // init & update, so could be calculated there.
    
    // Calculate relative availability correction, so calls from vectorUpdate,
    // etc., will have a mean of 1.0.
    double sumRelativeAvailability = 0.0;
    for (std::list<Host::Human>::iterator h = population.begin(); h != population.end(); ++h){
        sumRelativeAvailability += h->perHostTransmission.relativeAvailabilityAge (h->getAgeInYears());
    }
    ageCorrectionFactor = populationSize / sumRelativeAvailability;     // 1 / mean-rel-avail
    if( sumRelativeAvailability == 0.0 )
        // value should be unimportant when no humans are available, though inf/nan is not acceptable
        ageCorrectionFactor = 1.0;
}

void TransmissionModel::updateKappa (const std::list<Host::Human>& population) {
  // We calculate kappa for output and non-vector model, and kappaByAge for
  // the shared graphics.
  
  double sumWt_kappa= 0.0;
  double sumWeight  = 0.0;
  kappaByAge.assign (noOfAgeGroupsSharedMem, 0.0);
  nByAge.assign (noOfAgeGroupsSharedMem, 0);
  numTransmittingHumans = 0;
  
  for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
    double t = h->perHostTransmission.relativeAvailabilityHetAge(h->getAgeInYears());
    sumWeight += t;
    t *= h->probTransmissionToMosquito();
    sumWt_kappa += t;
    if( t > 0.0 )
        ++numTransmittingHumans;
    
    // kappaByAge and nByAge are used in the screensaver only
    Monitoring::AgeGroup ag = h->ageGroup();
    kappaByAge[ag.i()] += t;
    ++nByAge[ag.i()];
  }
  
  
  int tmod = TimeStep::simulation % TimeStep::stepsPerYear;
  int t1mod = (TimeStep::simulation-TimeStep(1)) % TimeStep::stepsPerYear;
  if( population.empty() ){     // this is valid
      kappa[tmod] = 0.0;        // no humans: no infectiousness
  } else {
    if ( !(sumWeight > DBL_MIN * 10.0) ){       // if approx. eq. 0, negative or an NaN
        ostringstream msg;
        msg<<"sumWeight is invalid: "<<sumWeight<<", "<<sumWt_kappa<<", "<<population.size();
        throw runtime_error(msg.str());
    }
    kappa[tmod] = sumWt_kappa / sumWeight;
  }
  
  //Calculate time-weighted average of kappa
  _sumAnnualKappa += kappa[tmod] * initialisationEIR[t1mod];
  if (tmod == 0) {
      // if annualEIR == 0.0 (or an NaN), we just get some nonsense output like inf or nan.
      // This is a better solution than printing a warning no-one will see and outputting 0.
      _annualAverageKappa = _sumAnnualKappa / annualEIR;
      _sumAnnualKappa = 0.0;
  }
  
  // Shared graphics: report infectiousness
  if (TimeStep::simulation % 6 ==  0) {
    for (size_t i = 0; i < noOfAgeGroupsSharedMem; i++)
      kappaByAge[i] /= nByAge[i];
    util::SharedGraphics::copyKappa(&kappaByAge[0]);
  }
  
  // Sum up innoculations this timestep
  double timeStepTotal = 0.0;
  for (size_t group = 0; group < timeStepEntoInnocs.size(); ++group) {
    timeStepTotal += timeStepEntoInnocs[group];
    innoculationsPerAgeGroup[group] += timeStepEntoInnocs[group];
    // Reset to zero:
    timeStepEntoInnocs[group] = 0.0;
  }
  innoculationsPerDayOfYear[tmod] = timeStepTotal / timeStepNumEntoInnocs;

  BSSInitialisationEIR += initialisationEIR[tmod];
  BSSInnoculationsPerDayOfYear +=innoculationsPerDayOfYear[tmod];
  BSSTimesteps++;

  timeStepNumEntoInnocs = 0;
}

double TransmissionModel::getEIR (OM::Transmission::PerHostTransmission& host, double ageYears, OM::Monitoring::AgeGroup ageGroup) {
  /* For the NonVector model, the EIR should just be multiplied by the
   * availability. For the Vector model, the availability is also required
   * for internal calculations, but again the EIR should be multiplied by the
   * availability. */
  double EIR = calculateEIR (host, ageYears);
  
  timeStepEntoInnocs[ageGroup.i()] += EIR;
  timeStepNumEntoInnocs ++;
  util::streamValidate( EIR );
  return EIR;
}

void TransmissionModel::summarize (Monitoring::Survey& survey) {
  survey.setNumTransmittingHosts(kappa[TimeStep::simulation % TimeStep::stepsPerYear]);
  survey.setAnnualAverageKappa(_annualAverageKappa);
  
  survey.setInnoculationsPerAgeGroup (innoculationsPerAgeGroup);        // Array contents must be copied.
  innoculationsPerAgeGroup.assign (innoculationsPerAgeGroup.size(), 0.0);

  survey.set_Vector_EIR_Input (BSSInitialisationEIR/(double)BSSTimesteps);
  survey.set_Vector_EIR_Simulated (BSSInnoculationsPerDayOfYear/(double)BSSTimesteps);

  BSSInitialisationEIR = 0.0;
  BSSInnoculationsPerDayOfYear = 0.0;
  BSSTimesteps = 0;
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
    BSSInitialisationEIR & stream;
    BSSInnoculationsPerDayOfYear & stream;
    BSSTimesteps & stream;
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
    BSSInitialisationEIR & stream;
    BSSInnoculationsPerDayOfYear & stream;
    BSSTimesteps & stream;
}

} }
