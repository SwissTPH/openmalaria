/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2013 Swiss Tropical and Public Health Institute 
 * Copyright (C) 2005-2013 Liverpool School Of Tropical Medicine
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

#include "Transmission/TransmissionModel.h"
#include "Monitoring/Surveys.h"
#include "PopulationStats.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "Population.h"
#include "interventions/Interventions.h"
#include <schema/scenario.h>

#include <string>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <gsl/gsl_cdf.h>

namespace OM { namespace Host {
    using namespace OM::util;
    int Human::_ylagLen = 0;

// -----  Static functions  -----

void Human::initHumanParameters( const Parameters& parameters, const scnXml::Scenario& scenario ) {    // static
    const scnXml::Model& model = scenario.getModel();
    // Init models used by humans:
    Transmission::PerHost::init( model.getHuman().getAvailabilityToMosquitoes() );
    InfectionIncidenceModel::init( parameters );
    WithinHost::WithinHostModel::init( parameters, scenario );
    Clinical::ClinicalModel::init( parameters, model, scenario.getHealthSystem() );
    _ylagLen = TimeStep::intervalsPer5Days.asInt() * 4;
}

void Human::clear() {   // static clear
  Clinical::ClinicalModel::cleanup();
  Transmission::PerHost::cleanup();
}


// -----  Non-static functions: creation/destruction, checkpointing  -----

// Create new human
Human::Human(Transmission::TransmissionModel& tm, TimeStep dateOfBirth) :
    perHostTransmission(tm),
    withinHostModel(WithinHost::WithinHostModel::createWithinHostModel()),
    infIncidence(InfectionIncidenceModel::createModel()),
    _dateOfBirth(dateOfBirth),
    nextCtsDist(0),
    _probTransmissionToMosquito(0.0)
{
  // Initial humans are created at time 0 and may have DOB in past. Otherwise DOB must be now.
  assert( _dateOfBirth == TimeStep::simulation ||
      (TimeStep::simulation == TimeStep(0) && _dateOfBirth < TimeStep::simulation));
  
  _ylag.assign (_ylagLen, 0.0);
  
  
  /* Human heterogeneity; affects:
   * _comorbidityFactor (stored in PathogenesisModel)
   * _treatmentSeekingFactor (stored in CaseManagementModel)
   * availabilityFactor (stored in Transmission::PerHost)
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

bool Human::update(Transmission::TransmissionModel* transmissionModel, bool doUpdate) {
#ifdef WITHOUT_BOINC
    ++PopulationStats::humanUpdateCalls;
    if( doUpdate )
        ++PopulationStats::humanUpdates;
#endif
    TimeStep ageTimeSteps = TimeStep::simulation-_dateOfBirth;
    if (clinicalModel->isDead(ageTimeSteps))
        return true;
    
    if (doUpdate){
        util::streamValidate( ageTimeSteps.asInt() );
        double ageYears = ageTimeSteps.inYears();
        monitoringAgeGroup.update( ageYears );
        
        updateInfection(transmissionModel, ageYears);
        clinicalModel->update (*this, ageYears, ageTimeSteps);
        clinicalModel->updateInfantDeaths (ageTimeSteps);
    }
    return false;
}

void Human::addInfection(){
    withinHostModel->importInfection();
}

void Human::updateInfection(Transmission::TransmissionModel* transmissionModel, double ageYears){
    // Cache total density for infectiousness calculations
    _ylag[mod_nn(TimeStep::simulation.asInt(),_ylagLen)] = withinHostModel->getTotalDensity();
    
    double EIR = transmissionModel->getEIR( perHostTransmission, ageYears, monitoringAgeGroup );
    int nNewInfs = infIncidence->numNewInfections( *this, EIR );
    
    withinHostModel->update(nNewInfs, ageYears, _vaccine.getEfficacy(Vaccine::BSV));
}

bool Human::needsRedeployment( size_t effect_index, TimeStep maxAge ){
    map<size_t,TimeStep>::const_iterator it = lastDeployments.find( effect_index );
    if( it == lastDeployments.end() ){
        return true;  // no previous deployment
    }else{
        return it->second + maxAge <= TimeStep::simulation;
    }
}

void Human::deployVaccine( interventions::Deployment::Method method, Vaccine::Types type ){
    if( method == interventions::Deployment::TIMED ){
        _vaccine.vaccinate( type );
        if( Vaccine::reportFor( type ) )
            Monitoring::Surveys.getSurvey(isInAnyCohort()).reportMassVaccinations (getMonitoringAgeGroup(), 1);
    }else if( method == interventions::Deployment::CTS ){
        if ( _vaccine.getsEPIVaccination( type, TimeStep::simulation - _dateOfBirth ) ){
            _vaccine.vaccinate( type );
            if( Vaccine::reportFor( type ) )
                Monitoring::Surveys.getSurvey(isInAnyCohort()).reportEPIVaccinations (getMonitoringAgeGroup(), 1);
        }
    }else throw SWITCH_DEFAULT_EXCEPTION;
}

void Human::massDrugAdministration () {
    clinicalModel->massDrugAdministration (*this);
}

void Human::reportDeployment( interventions::Effect::Type type, interventions::Deployment::Method method ) const{
    if( method == interventions::Deployment::TIMED ){
        switch( type ){
            case interventions::Effect::GVI:
                Monitoring::Surveys.getSurvey(isInAnyCohort()).reportMassGVI( getMonitoringAgeGroup(), 1 );
                break;
            default:
                throw SWITCH_DEFAULT_EXCEPTION;
        }
    }else if( method == interventions::Deployment::CTS ){
        switch( type ){
            case interventions::Effect::GVI:
                //TODO(monitoring): report
                break;
            default:
                throw SWITCH_DEFAULT_EXCEPTION;
        }
    }else throw SWITCH_DEFAULT_EXCEPTION;
}

double Human::getAgeInYears() const{
    return (TimeStep::simulation - _dateOfBirth).inYears();
}


void Human::summarize() {
    // 5-day only, compatibility option:
    if( util::ModelOptions::option( util::REPORT_ONLY_AT_RISK ) &&
        clinicalModel->notAtRisk() ){
        // This modifies the denominator to treat the 4*5 day intervals
        // after a bout as 'not at risk' to match the IPTi trials
        return;
    }
    
    Monitoring::Survey& survey( Monitoring::Surveys.getSurvey( isInAnyCohort() ) );
    survey.reportHosts (getMonitoringAgeGroup(), 1);
    bool patent = withinHostModel->summarize (survey, getMonitoringAgeGroup());
    infIncidence->summarize (survey, getMonitoringAgeGroup());
    clinicalModel->summarize (survey, getMonitoringAgeGroup());
    
    if( patent ){
        removeFromCohorts( interventions::CohortSelectionEffect::REMOVE_AT_FIRST_INFECTION );
    }
}

void Human::addToCohort (size_t index){
    if( cohorts.count(index) > 0 ) return;	// nothing to do
    // Data accumulated between reports should be flushed. Currently all this
    // data remembers which survey it should go to or is reported immediately,
    // although episode reports still need to be flushed.
    flushReports();
    cohorts.insert(index);
    //TODO(monitoring): reporting is inappropriate
    Monitoring::Surveys.current->reportAddedToCohort( getMonitoringAgeGroup(), 1 );
}
void Human::removeFromCohort( size_t index ){
    if( cohorts.count(index) > 0 ){
        // Data should be flushed as with addToCohort().
        flushReports();
        cohorts.erase( index );
        //TODO(monitoring): reporting
        Monitoring::Surveys.current->reportRemovedFromCohort( getMonitoringAgeGroup(), 1 );
    }
}
void Human::removeFromCohorts( interventions::CohortSelectionEffect::RemoveAtCode code ){
    const vector<size_t>& removeAtList = interventions::CohortSelectionEffect::removeAtIds[code];
    for( vector<size_t>::const_iterator it = removeAtList.begin(), end = removeAtList.end(); it != end; ++it ){
        removeFromCohort( *it );    // only does anything if in cohort
    }
}

void Human::flushReports (){
    clinicalModel->flushReports();
}

void Human::updateInfectiousness() {
  /* This model (often referred to as the gametocyte model) was designed for
  5-day timesteps. We use the same model (sampling 10, 15 and 20 days ago)
  for 1-day timesteps to avoid having to design and analyse a new model.
  Description: AJTMH pp.32-33 */
  TimeStep ageTimeSteps=TimeStep::simulation-_dateOfBirth;
  if (ageTimeSteps.inDays() <= 20 || TimeStep::simulation.inDays() <= 20){
    // We need at least 20 days history (_ylag) to calculate infectiousness;
    // assume no infectiousness if we don't have this history.
    // Note: human not updated on DOB so age must be >20 days.
    return;
  }
  
  //Infectiousness parameters: see AJTMH p.33, tau=1/sigmag**2 
  static const double beta1=1.0;
  static const double beta2=0.46;
  static const double beta3=0.17;
  static const double tau= 0.066;
  static const double mu= -8.1;
  
  // Take weighted sum of total asexual blood stage density 10, 15 and 20 days
  // before. We have 20 days history, so use mod_nn:
  int firstIndex = TimeStep::simulation.asInt()-2*TimeStep::intervalsPer5Days.asInt() + 1;
  double x = beta1 * _ylag[mod_nn(firstIndex, _ylagLen)]
           + beta2 * _ylag[mod_nn(firstIndex-TimeStep::intervalsPer5Days.asInt(), _ylagLen)]
           + beta3 * _ylag[mod_nn(firstIndex-2*TimeStep::intervalsPer5Days.asInt(), _ylagLen)];
  if (x < 0.001){
    _probTransmissionToMosquito = 0.0;
    return;
  }
  
  double zval=(log(x)+mu)/sqrt(1.0/tau);
  double pone = gsl_cdf_ugaussian_P(zval);
  double transmit=(pone*pone);
  //transmit has to be between 0 and 1
  transmit=std::max(transmit, 0.0);
  transmit=std::min(transmit, 1.0);
  
  //    Include here the effect of transmission-blocking vaccination
  _probTransmissionToMosquito = transmit*(1.0-_vaccine.getEfficacy( Vaccine::TBV ));
  util::streamValidate( _probTransmissionToMosquito );
}

} }
