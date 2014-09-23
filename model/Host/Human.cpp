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
#include "interventions/InterventionManager.hpp"
#include "schema/scenario.h"

namespace OM { namespace Host {
    using namespace OM::util;
    using namespace Monitoring;
    using interventions::ComponentId;
    
    bool opt_trans_het = false, opt_comorb_het = false, opt_treat_het = false,
            opt_trans_treat_het = false, opt_comorb_treat_het = false,
            opt_comorb_trans_het = false, opt_triple_het = false,
            opt_report_only_at_risk = false;
    
    //HACK(drug mon)
    ofstream monDrug, monFake;
    ComponentId drugMonId(0 /*lazy initialisation*/);

// -----  Static functions  -----

void Human::init( const Parameters& parameters, const scnXml::Scenario& scenario ){    // static
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
    Clinical::ClinicalModel::init( parameters, scenario );
}

void Human::init2( const scnXml::Monitoring& monitoring ){
    if( monitoring.getDrugConcentration().present() ){
        monDrug.open( monitoring.getDrugConcentration().get().getFile().c_str(), ios::out );
        drugMonId = interventions::InterventionManager::getComponentId( monitoring.getDrugConcentration().get().getCohort() );
    }
}


// -----  Non-static functions: creation/destruction, checkpointing  -----

// Create new human
Human::Human(Transmission::TransmissionModel& tm, SimTime dateOfBirth) :
    infIncidence(InfectionIncidenceModel::createModel()),
    m_DOB(dateOfBirth),
    m_cohortSet(0),
    nextCtsDist(0)
{
  // Initial humans are created at time 0 and may have DOB in past. Otherwise DOB must be now.
  assert( m_DOB == sim::now1() || (sim::now1() == sim::zero() && m_DOB < sim::now1()) );
  
  
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

Human::Human(SimTime dateOfBirth) :
    withinHostModel(0),
    infIncidence(0),
    clinicalModel(0),
    m_DOB(dateOfBirth),
    m_cohortSet(0),
    nextCtsDist(0)
{}

void Human::destroy() {
  delete infIncidence;
  delete withinHostModel;
  delete clinicalModel;
}


// -----  Non-static functions: per-time-step update  -----

bool Human::update(Transmission::TransmissionModel* transmissionModel, bool doUpdate) {
#ifdef WITHOUT_BOINC
    ++PopulationStats::humanUpdateCalls;
    if( doUpdate )
        ++PopulationStats::humanUpdates;
#endif
    SimTime age = getAge();
    if (clinicalModel->isDead(age))
        return true;
    
    if (doUpdate){
        util::streamValidate( age.raw() );
        double ageYears = age.inYears();
        monitoringAgeGroup.update( ageYears );
        // check sub-pop expiry
        for( map<ComponentId,SimTime>::iterator expIt =
            m_subPopExp.begin(), expEnd = m_subPopExp.end(); expIt != expEnd; )
        {
            //  We use >= to test membership, < is inverse (see comment on m_subPopExp)
            if( expIt->second < sim::now1() ){
                // don't flush reports
                // report removal due to expiry
                Survey::current().addInt(Report::MI_N_SP_REM_TOO_OLD, *this, 1 );
                m_cohortSet = Survey::updateCohortSet( m_cohortSet, expIt->first, false );
                // erase element, but continue iteration (note: this is simpler in C++11)
                map<ComponentId,SimTime>::iterator toErase = expIt;
                ++expIt;
                m_subPopExp.erase( toErase );
            }else{
                ++expIt;
            }
        }
        double EIR = transmissionModel->getEIR( *this, ageYears, monitoringAgeGroup );
        int nNewInfs = infIncidence->numNewInfections( *this, EIR );
        
        ofstream& mon = isInSubPop(drugMonId) ? monDrug : monFake;
        withinHostModel->update(nNewInfs, ageYears, _vaccine.getFactor(interventions::Vaccine::BSV), mon);
        
        clinicalModel->update( *this, ageYears, age == sim::oneTS() );
        clinicalModel->updateInfantDeaths( age );
    }
    return false;
}

void Human::addInfection(){
    withinHostModel->importInfection();
}

void Human::clearImmunity(){
    withinHostModel->clearImmunity();
}


void Human::summarize() {
    // 5-day only, compatibility option:
    if( opt_report_only_at_risk && clinicalModel->notAtRisk() ){
        // This modifies the denominator to treat the 4*5 day intervals
        // after a bout as 'not at risk' to match the IPTi trials
        return;
    }
    
    Survey::current().addInt( Report::MI_HOSTS, *this, 1)
        .addDouble( Report::MD_AGE, *this, getAgeInYears() );
    bool patent = withinHostModel->summarize (*this);
    infIncidence->summarize (*this);
    
    if( patent ){
        // this should happen after all other reporting!
        removeFirstEvent( interventions::SubPopRemove::ON_FIRST_INFECTION );
    }
}

void Human::reportDeployment( ComponentId id, SimTime duration ){
    if( duration <= sim::zero() ) return; // nothing to do
    m_subPopExp[id] = sim::now1() + duration;
    m_cohortSet = Survey::updateCohortSet( m_cohortSet, id, true );
}
void Human::removeFirstEvent( interventions::SubPopRemove::RemoveAtCode code ){
    const vector<ComponentId>& removeAtList = interventions::removeAtIds[code];
    for( vector<ComponentId>::const_iterator it = removeAtList.begin(), end = removeAtList.end(); it != end; ++it ){
        SubPopT::iterator expIt = m_subPopExp.find( *it );
        if( expIt != m_subPopExp.end() ){
            //  We use >= to test membership (see comment on m_subPopExp)
            if( expIt->second >= sim::now1() ){
                // removeFirstEvent() is used for onFirstBout, onFirstTreatment
                // and onFirstInfection cohort options. Health system memory must
                // be reset for this to work properly; in theory the memory should
                // be independent for each cohort, but this is a usable approximation.
                flushReports();     // reset HS memory
                
                // report removal due to first infection/bout/treatment
                Survey::current().addInt(Report::MI_N_SP_REM_FIRST_EVENT, *this, 1 );
            }
            m_cohortSet = Survey::updateCohortSet( m_cohortSet, expIt->first, false );
            // remove (affects reporting, restrictToSubPop and cumulative deployment):
            m_subPopExp.erase( expIt );
        }
    }
}

void Human::flushReports (){
    clinicalModel->flushReports();
}

} }
