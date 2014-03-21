/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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
#include "WithinHost/WHInterface.h"

#include "Transmission/TransmissionModel.h"
#include "Monitoring/Survey.h"
#include "PopulationStats.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "Population.h"
#include "interventions/Cohort.h"
#include <schema/scenario.h>

namespace OM { namespace Host {
    using namespace OM::util;
    using namespace Monitoring;
    
    bool opt_trans_het = false, opt_comorb_het = false, opt_treat_het = false,
            opt_trans_treat_het = false, opt_comorb_treat_het = false,
            opt_comorb_trans_het = false, opt_triple_het = false,
            opt_report_only_at_risk = false;

// -----  Static functions  -----

void Human::initHumanParameters( const Parameters& parameters, const scnXml::Scenario& scenario ) {    // static
    opt_trans_het = util::ModelOptions::option (util::TRANS_HET);
    opt_comorb_het = util::ModelOptions::option (util::COMORB_HET);
    opt_treat_het = util::ModelOptions::option (util::TREAT_HET);
    opt_trans_treat_het = util::ModelOptions::option (util::TRANS_TREAT_HET);
    opt_comorb_treat_het = util::ModelOptions::option (util::COMORB_TREAT_HET);
    opt_comorb_trans_het = util::ModelOptions::option (util::COMORB_TRANS_HET);
    opt_triple_het = util::ModelOptions::option (util::TRIPLE_HET);
    opt_report_only_at_risk = util::ModelOptions::option( util::REPORT_ONLY_AT_RISK );
    
    const scnXml::Model& model = scenario.getModel();
    // Init models used by humans:
    Transmission::PerHost::init( model.getHuman().getAvailabilityToMosquitoes() );
    InfectionIncidenceModel::init( parameters );
    WithinHost::WHInterface::init( parameters, scenario );
    Clinical::ClinicalModel::init( parameters, model, scenario.getHealthSystem() );
}

void Human::clear() {   // static clear
  Clinical::ClinicalModel::cleanup();
  Transmission::PerHost::cleanup();
}


// -----  Non-static functions: creation/destruction, checkpointing  -----

// Create new human
Human::Human(Transmission::TransmissionModel& tm, TimeStep dateOfBirth) :
    perHostTransmission(tm),
    infIncidence(InfectionIncidenceModel::createModel()),
    _dateOfBirth(dateOfBirth),
    nextCtsDist(0)
{
  // Initial humans are created at time 0 and may have DOB in past. Otherwise DOB must be now.
  assert( _dateOfBirth == TimeStep::simulation ||
      (TimeStep::simulation == TimeStep(0) && _dateOfBirth < TimeStep::simulation));
  
  
  /* Human heterogeneity; affects:
   * _comorbidityFactor (stored in PathogenesisModel)
   * _treatmentSeekingFactor (stored in CaseManagementModel)
   * availabilityFactor (stored in Transmission::PerHost)
   */
  double comorbidityFactor = 1.0;
  double treatmentSeekingFactor = 1.0;
  double availabilityFactor = 1.0;
  
  if (opt_trans_het) {
    availabilityFactor=0.2;
    if (random::uniform_01() < 0.5) {
      availabilityFactor=1.8;
    }
  }
  if (opt_comorb_het) {
    comorbidityFactor=0.2;
    if (random::uniform_01() < 0.5) {
      comorbidityFactor=1.8;
    }   
  }
  if (opt_treat_het) {
    treatmentSeekingFactor=0.2;
    if (random::uniform_01() < 0.5) {            
      treatmentSeekingFactor=1.8;
    }   
  }
  if (opt_trans_treat_het) {
    treatmentSeekingFactor=0.2;
    availabilityFactor=1.8;
    if (random::uniform_01()<0.5) {
      treatmentSeekingFactor=1.8;
      availabilityFactor=0.2;
    }
  } else if (opt_comorb_treat_het) {
    if (random::uniform_01()<0.5) {
      comorbidityFactor=1.8;
      treatmentSeekingFactor=0.2;
    } else {
      comorbidityFactor=0.2;
      treatmentSeekingFactor=1.8;
    }
  } else if (opt_comorb_trans_het) {
    availabilityFactor=1.8;
    comorbidityFactor=1.8;
    if (random::uniform_01()<0.5) {
      availabilityFactor=0.2;
      comorbidityFactor=0.2;
    }
  } else if (opt_triple_het) {
    availabilityFactor=1.8;
    comorbidityFactor=1.8;
    treatmentSeekingFactor=0.2;
    if (random::uniform_01()<0.5) {
      availabilityFactor=0.2;
      comorbidityFactor=0.2;
      treatmentSeekingFactor=1.8;
    }
  }
  withinHostModel = WithinHost::WHInterface::createWithinHostModel( comorbidityFactor );
  perHostTransmission.initialise (tm, availabilityFactor * infIncidence->getAvailabilityFactor(1.0));
  clinicalModel = Clinical::ClinicalModel::createClinicalModel (treatmentSeekingFactor);
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
        
        double EIR = transmissionModel->getEIR( perHostTransmission, ageYears, monitoringAgeGroup );
        int nNewInfs = infIncidence->numNewInfections( *this, EIR );
        
        withinHostModel->update(nNewInfs, ageYears, _vaccine.getFactor(interventions::Vaccine::BSV));
        
        clinicalModel->update (*this, ageYears, ageTimeSteps);
        clinicalModel->updateInfantDeaths (ageTimeSteps);
    }
    return false;
}

void Human::addInfection(){
    withinHostModel->importInfection();
}

void Human::clearImmunity(){
    withinHostModel->clearImmunity();
}

bool Human::needsRedeployment( interventions::ComponentId cumCovId, TimeStep maxAge ){
    map<interventions::ComponentId,TimeStep>::const_iterator it = lastDeployments.find( cumCovId );
    if( it == lastDeployments.end() ){
        return true;  // no previous deployment
    }else{
        return it->second + maxAge <= TimeStep::simulation;
    }
}


void Human::summarize() {
    // 5-day only, compatibility option:
    if( opt_report_only_at_risk && clinicalModel->notAtRisk() ){
        // This modifies the denominator to treat the 4*5 day intervals
        // after a bout as 'not at risk' to match the IPTi trials
        return;
    }
    
    Survey::current().addInt( Report::MI_HOSTS, *this, 1);
    bool patent = withinHostModel->summarize (*this);
    infIncidence->summarize (*this);
    
    if( patent ){
        // this should happen after all other reporting!
        removeFromCohorts( interventions::Cohort::REMOVE_AT_FIRST_INFECTION );
    }
}

void Human::addToCohort (interventions::ComponentId id, ReportMeasureI reportMeasure){
    if( cohorts.count(id) > 0 ) return;	// nothing to do
    
    // Data accumulated between reports should be flushed. Currently all this
    // data remembers which survey it should go to or is reported immediately,
    // although episode reports still need to be flushed.
    flushReports();
    cohorts.insert(id);
    //TODO(monitoring): reporting is inappropriate
    Survey::current().addInt( reportMeasure, *this, 1 );
}
void Human::removeFromCohort( interventions::ComponentId id ){
    if( cohorts.count(id) <= 0 ) return;        // nothing to do
    
    // Data should be flushed as with addToCohort().
    flushReports();
    cohorts.erase( id );
    //TODO(monitoring): reporting
    Survey::current().addInt(Report::MI_NUM_REMOVED_COHORT, *this, 1 );
}
void Human::removeFromCohorts( interventions::Cohort::RemoveAtCode code ){
    const vector<interventions::ComponentId>& removeAtList = interventions::CohortSelectionComponent::removeAtIds[code];
    for( vector<interventions::ComponentId>::const_iterator it = removeAtList.begin(), end = removeAtList.end(); it != end; ++it ){
        removeFromCohort( *it );    // only does anything if in cohort
    }
}

void Human::flushReports (){
    clinicalModel->flushReports();
}

} }
