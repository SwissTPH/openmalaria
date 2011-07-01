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

#include "Host/Human.h"

#include "Host/InfectionIncidenceModel.h"
#include "Clinical/ClinicalModel.h"
#include "WithinHost/DescriptiveIPTWithinHost.h"        // only for summarizing

#include "inputData.h"
#include "Transmission/TransmissionModel.h"
#include "Monitoring/Surveys.h"
#include "PopulationStats.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "Interventions.h"

#include <string>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <gsl/gsl_cdf.h>

namespace OM { namespace Host {
    using namespace OM::util;
    int Human::_ylagLen = 0;
    bool Human::cohortFirstBoutOnly = false;
    bool Human::cohortFirstTreatmentOnly = false;
    bool Human::cohortFirstInfectionOnly = false;

// -----  Static functions  -----

void Human::initHumanParameters () {    // static
    // Init models used by humans:
    Transmission::PerHostTransmission::init();
    InfectionIncidenceModel::init();
    WithinHost::WithinHostModel::init();
    Clinical::ClinicalModel::init();
    _ylagLen = TimeStep::intervalsPer5Days.asInt() * 4;
    
    cohortFirstBoutOnly = InputData().getMonitoring().getFirstBoutOnly();
    cohortFirstTreatmentOnly = InputData().getMonitoring().getFirstTreatmentOnly();
    cohortFirstInfectionOnly = InputData().getMonitoring().getFirstInfectionOnly();
}

void Human::clear() {   // static clear
  Clinical::ClinicalModel::cleanup();
  Transmission::PerHostTransmission::cleanup();
}


// -----  Non-static functions: creation/destruction, checkpointing  -----

// Create new human
Human::Human(Transmission::TransmissionModel& tm, TimeStep dateOfBirth) :
    perHostTransmission(),
    withinHostModel(WithinHost::WithinHostModel::createWithinHostModel()),
    infIncidence(InfectionIncidenceModel::createModel()),
    _dateOfBirth(dateOfBirth),
    nextCtsDist(0),
    _inCohort(false),
    _probTransmissionToMosquito(0.0)
{
  if( _dateOfBirth != TimeStep::simulation1() &&
      (TimeStep::simulation != TimeStep(0) || _dateOfBirth > TimeStep::simulation))
  {
    // Initial humans are created at time 0 and may have DOB in past. Otherwise DOB must be now.
    throw util::traced_exception ("Invalid date of birth!");
  }
  
  _ylag.assign (_ylagLen, 0.0);
  
  
  /* Human heterogeneity; affects:
   * _comorbidityFactor (stored in PathogenesisModel)
   * _treatmentSeekingFactor (stored in CaseManagementModel)
   * availabilityFactor (stored in Transmission::PerHostTransmission)
   */
  double _comorbidityFactor = 1.0;
  double _treatmentSeekingFactor = 1.0;
  double availabilityFactor = 1.0;
  
  if (util::ModelOptions::option (util::TRANS_HET)) {
    availabilityFactor=0.2;
    if (random::uniform_01() < 0.5) {
      availabilityFactor=1.8;
    }
  }
  if (util::ModelOptions::option (util::COMORB_HET)) {
    _comorbidityFactor=0.2;
    if (random::uniform_01() < 0.5) {
      _comorbidityFactor=1.8;
    }   
  }
  if (util::ModelOptions::option (util::TREAT_HET)) {
    _treatmentSeekingFactor=0.2;
    if (random::uniform_01() < 0.5) {            
      _treatmentSeekingFactor=1.8;
    }   
  }
  if (util::ModelOptions::option (util::TRANS_TREAT_HET)) {
    _treatmentSeekingFactor=0.2;
    availabilityFactor=1.8;
    if (random::uniform_01()<0.5) {
      _treatmentSeekingFactor=1.8;
      availabilityFactor=0.2;
    }
  } else if (util::ModelOptions::option (util::COMORB_TREAT_HET)) {
    if (random::uniform_01()<0.5) {
      _comorbidityFactor=1.8;
      _treatmentSeekingFactor=0.2;
    } else {
      _comorbidityFactor=0.2;
      _treatmentSeekingFactor=1.8;
    }
  } else if (util::ModelOptions::option (util::COMORB_TRANS_HET)) {
    availabilityFactor=1.8;
    _comorbidityFactor=1.8;
    if (random::uniform_01()<0.5) {
      availabilityFactor=0.2;
      _comorbidityFactor=0.2;
    }
  } else if (util::ModelOptions::option (util::TRIPLE_HET)) {
    availabilityFactor=1.8;
    _comorbidityFactor=1.8;
    _treatmentSeekingFactor=0.2;
    if (random::uniform_01()<0.5) {
      availabilityFactor=0.2;
      _comorbidityFactor=0.2;
      _treatmentSeekingFactor=1.8;
    }
  }
  perHostTransmission.initialise (tm, availabilityFactor * infIncidence->getAvailabilityFactor(1.0));
  clinicalModel = Clinical::ClinicalModel::createClinicalModel (_comorbidityFactor, _treatmentSeekingFactor);
}

void Human::destroy() {
  delete infIncidence;
  delete withinHostModel;
  delete clinicalModel;
}


// -----  Non-static functions: per-timestep update  -----

bool Human::update(const OM::Population& population, Transmission::TransmissionModel* transmissionModel, bool doUpdate) {
#ifdef WITHOUT_BOINC
    ++PopulationStats::humanUpdateCalls;
    if( doUpdate )
        ++PopulationStats::humanUpdates;
#endif
    TimeStep ageTimeSteps = TimeStep::simulation1()-_dateOfBirth;
    if (clinicalModel->isDead(ageTimeSteps))
        return true;
    
    if (doUpdate){
        util::streamValidate( ageTimeSteps.asInt() );
        double ageYears = ageTimeSteps.inYears();
        monitoringAgeGroup.update( ageYears );
        
        updateInterventionStatus(population);
        updateInfection(transmissionModel, ageYears);
        clinicalModel->update (*this, ageYears, ageTimeSteps);
        clinicalModel->updateInfantDeaths (ageTimeSteps);
        // set _probTransmissionToMosquito to next time-step's value
        _probTransmissionToMosquito = calcProbTransmissionToMosquito ();
    }
    return false;
}

void Human::addInfection(){
    withinHostModel->newInfection(TimeStep::simulation);
}

void Human::updateInfection(Transmission::TransmissionModel* transmissionModel, double ageYears){
    double EIR = transmissionModel->getEIR( perHostTransmission, ageYears, monitoringAgeGroup );
    int numInf = infIncidence->numNewInfections( *this, EIR );
    for (int i=1;i<=numInf; ++i) {
        withinHostModel->newInfection(TimeStep::simulation);
    }
    
    // Cache total density for infectiousness calculations
    _ylag[TimeStep::simulation1().asInt()%_ylagLen]=withinHostModel->getTotalDensity();
    
    withinHostModel->calculateDensities(ageYears, _vaccine.getBSVEfficacy());
}

void Human::updateInterventionStatus(const OM::Population& population) {
    if (TimeStep::interventionPeriod >= TimeStep(0)) {
        TimeStep ageTimeSteps = TimeStep::simulation1()-_dateOfBirth;
        //HACK
        InterventionManager::getSingleton().deployCts(population, *this, ageTimeSteps, nextCtsDist);
    }
}


void Human::massVaccinate (const OM::Population&) {
    _vaccine.vaccinate(TimeStep::simulation);
    Monitoring::Surveys.getSurvey(_inCohort).reportMassVaccinations (getMonitoringAgeGroup(), 1);
}
void Human::ctsVaccinate (const OM::Population&) {
    if ( _vaccine.doCtsVaccination( TimeStep::simulation1() - _dateOfBirth ) ){
        _vaccine.vaccinate(TimeStep::simulation);
        Monitoring::Surveys.getSurvey(_inCohort).reportEPIVaccinations (getMonitoringAgeGroup(), 1);
    }
}

void Human::IPTiTreatment (const OM::Population&) {
  withinHostModel->IPTiTreatment (getMonitoringAgeGroup(), _inCohort);
}
void Human::deployIptDose (const OM::Population&) {
    withinHostModel->deployIptDose( getMonitoringAgeGroup(), _inCohort );
}

void Human::massDrugAdministration (const OM::Population&) {
    clinicalModel->massDrugAdministration (*this);
}

void Human::massITN (const OM::Population& population){
    perHostTransmission.setupITN (TimeStep::simulation, population.transmissionModel());
    Monitoring::Surveys.getSurvey(_inCohort).reportMassITNs( getMonitoringAgeGroup(), 1 );
}
void Human::ctsITN (const OM::Population& population){
    perHostTransmission.setupITN (TimeStep::simulation1(), population.transmissionModel());
    Monitoring::Surveys.getSurvey(_inCohort).reportEPI_ITNs( getMonitoringAgeGroup(), 1 );
}

void Human::massIRS (const OM::Population&) {
    perHostTransmission.setupIRS ();
    Monitoring::Surveys.getSurvey(_inCohort).reportMassIRS( getMonitoringAgeGroup(), 1 );
}

void Human::massVA (const OM::Population&) {
    perHostTransmission.setupVA ();
    Monitoring::Surveys.getSurvey(_inCohort).reportMassVA( getMonitoringAgeGroup(), 1 );
}

bool Human::hasVaccineProtection(TimeStep maxInterventionAge) const{
    return _vaccine.hasProtection(maxInterventionAge);
}
bool Human::hasIPTiProtection(TimeStep maxInterventionAge) const{
    return withinHostModel->hasIPTiProtection(maxInterventionAge);
}
bool Human::hasITNProtection(TimeStep maxInterventionAge) const{
    return perHostTransmission.getITN().timeOfDeployment() + maxInterventionAge > TimeStep::simulation;
}
bool Human::hasIRSProtection(TimeStep maxInterventionAge) const{
    return perHostTransmission.hasIRSProtection(maxInterventionAge);
}
bool Human::hasVAProtection(TimeStep maxInterventionAge) const{
    return perHostTransmission.hasVAProtection(maxInterventionAge);
}

double Human::getAgeInYears() const{
    return (TimeStep::simulation - _dateOfBirth).inYears();
}
double Human::getAgeInYears1() const{
    return (TimeStep::simulation1() - _dateOfBirth).inYears();
}


void Human::summarize() {
    // 5-day only, compatibility option:
    if( util::ModelOptions::option( util::REPORT_ONLY_AT_RISK ) &&
        clinicalModel->notAtRisk() ){
        // This modifies the denominator to treat the 4*5 day intervals
        // after a bout as 'not at risk' to match the IPTi trials
        return;
    }
    
    Monitoring::Survey& survey( Monitoring::Surveys.getSurvey( _inCohort ) );
    survey.reportHosts (getMonitoringAgeGroup(), 1);
    bool patent = withinHostModel->summarize (survey, getMonitoringAgeGroup());
    infIncidence->summarize (survey, getMonitoringAgeGroup());
    clinicalModel->summarize (survey, getMonitoringAgeGroup());
    
    if( cohortFirstInfectionOnly && patent ){
        removeFromCohort();
    }
}

void Human::addToCohort (const OM::Population&){
    if( _inCohort ) return;	// nothing to do
    // Data accumulated between reports should be flushed. Currently all this
    // data remembers which survey it should go to or is reported immediately,
    // although episode reports still need to be flushed.
    flushReports();
    _inCohort = true;
    Monitoring::Surveys.current->reportAddedToCohort( getMonitoringAgeGroup(), 1 );
}
void Human::removeFromCohort(){
    if( _inCohort ){
        // Data should be flushed as with addToCohort().
        flushReports();
        _inCohort = false;
        Monitoring::Surveys.current->reportRemovedFromCohort( getMonitoringAgeGroup(), 1 );
    }
}


void Human::flushReports (){
    clinicalModel->flushReports();
}


double Human::calcProbTransmissionToMosquito() const {
  /* This model (often referred to as the gametocyte model) was designed for
  5-day timesteps. We use the same model (sampling 10, 15 and 20 days ago)
  for 1-day timesteps to avoid having to design and analyse a new model.
  Description: AJTMH pp.32-33 */
  TimeStep ageTimeSteps=TimeStep::simulation1()-_dateOfBirth;
  if (ageTimeSteps.inDays() <= 20 || TimeStep::simulation1().inDays() <= 20)
    return 0.0;
  
  //Infectiousness parameters: see AJTMH p.33, tau=1/sigmag**2 
  static const double beta1=1.0;
  static const double beta2=0.46;
  static const double beta3=0.17;
  static const double tau= 0.066;
  static const double mu= -8.1;
  
  // Take weighted sum of total asexual blood stage density 10, 15 and 20 days before.
  // These values are one timestep more recent than that, however the calculated
  // value is not used until the next timestep when then ages would be correct.
  // Min TimeStep::simulation is 20 days, so LHS of '%' operator is never negative.
  int firstIndex = TimeStep::simulation1().asInt()-2*TimeStep::intervalsPer5Days.asInt()+1;
  double x = beta1 * _ylag[firstIndex % _ylagLen]
           + beta2 * _ylag[(firstIndex-TimeStep::intervalsPer5Days.asInt()) % _ylagLen]
           + beta3 * _ylag[(firstIndex-2*TimeStep::intervalsPer5Days.asInt()) % _ylagLen];
  if (x < 0.001)
    return 0.0;
  
  double zval=(log(x)+mu)/sqrt(1.0/tau);
  double pone = gsl_cdf_ugaussian_P(zval);
  double transmit=(pone*pone);
  //transmit has to be between 0 and 1
  transmit=std::max(transmit, 0.0);
  transmit=std::min(transmit, 1.0);
  
  //    Include here the effect of transmission-blocking vaccination
  double ret = transmit*(1.0-_vaccine.getTBVEfficacy());
  util::streamValidate( ret );
  return ret;
}

} }
