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

#include "Host/HumanHet.hpp"
#include "Host/InfectionIncidenceModel.h"
#include "Clinical/ClinicalModel.h"
#include "WithinHost/WHInterface.h"

#include "Transmission/TransmissionModel.h"
#include "PopulationStats.h"
#include "util/ModelOptions.h"
#include "util/vectors.h"
#include "util/StreamValidator.h"
#include "Population.h"
#include "interventions/InterventionManager.hpp"
#include "mon/reporting.h"
#include "schema/scenario.h"

namespace OM { namespace Host {
    using namespace OM::util;
    using Monitoring::Survey;
    using interventions::ComponentId;
    
    bool opt_report_only_at_risk = false;
    
    // Only required for a drug monitoring HACK and could be removed:
    ofstream monDrug, monFake;
    ComponentId drugMonId(0 /*lazy initialisation*/);

// -----  Static functions  -----

void Human::init( const Parameters& parameters, const scnXml::Scenario& scenario ){    // static
    HumanHet::init();
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
#ifdef WITHOUT_BOINC
        monDrug.open( monitoring.getDrugConcentration().get().getFile().c_str(), ios::out );
        drugMonId = interventions::InterventionManager::getComponentId( monitoring.getDrugConcentration().get().getCohort() );
#else
        // This feature is disabled in BOINC because it writes to a file and I
        // don't want another security feature needing review.
        std::cerr << "monitoring/drugConcentration: feature disabled in BOINC builds" << std::endl;
#endif
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
  assert( m_DOB == sim::nowOrTs1() || (sim::now() == sim::zero() && m_DOB < sim::now()) );
  
  HumanHet het = HumanHet::sample();
  withinHostModel = WithinHost::WHInterface::createWithinHostModel( het.comorbidityFactor );
  perHostTransmission.initialise (tm, het.availabilityFactor * infIncidence->getAvailabilityFactor(1.0));
  clinicalModel = Clinical::ClinicalModel::createClinicalModel (het.treatmentSeekingFactor);
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
    // For integer age checks we use age0 to e.g. get 73 steps comparing less than 1 year old
    SimTime age0 = age(sim::ts0());
    if (clinicalModel->isDead(age0))
        return true;
    
    if (doUpdate){
        util::streamValidate( age0.raw() );
        // Ages at  the end of the update period respectively. In most cases
        // the difference between this and age at the start is not especially
        // important in the model design, but since we parameterised with
        // ageYears1 we should stick with it.
        double ageYears1 = age(sim::ts1()).inYears();
        // monitoringAgeGroup is the group for the start of the time step.
        monitoringAgeGroup.update( age0 );
        // check sub-pop expiry
        for( map<ComponentId,SimTime>::iterator expIt =
            m_subPopExp.begin(), expEnd = m_subPopExp.end(); expIt != expEnd; )
        {
            if( !(expIt->second >= sim::ts0()) ){       // membership expired
                // don't flush reports
                // report removal due to expiry
                mon::reportMHI( mon::MHR_SUB_POP_REM_TOO_OLD, *this, 1 );
                m_cohortSet = Survey::updateCohortSet( m_cohortSet, expIt->first, false );
                // erase element, but continue iteration (note: this is simpler in C++11)
                map<ComponentId,SimTime>::iterator toErase = expIt;
                ++expIt;
                m_subPopExp.erase( toErase );
            }else{
                ++expIt;
            }
        }
        vector<double> EIR_per_genotype;        //TODO: avoid reallocating on every use!
        // ageYears1 used only in PerHost::relativeAvailabilityAge(); difference to age0 should be minor
        double EIR = transmissionModel->getEIR( *this, age0, ageYears1,
                monitoringAgeGroup, EIR_per_genotype );
        int nNewInfs = infIncidence->numNewInfections( *this, EIR );
        
        ofstream& mon = isInSubPop(drugMonId) ? monDrug : monFake;
        // ageYears1 used when medicating drugs (small effect) and in immunity model (which was parameterised for it)
        withinHostModel->update(nNewInfs, EIR_per_genotype, ageYears1,
                _vaccine.getFactor(interventions::Vaccine::BSV), mon);
        
        // ageYears1 used to get case fatality and sequelae probabilities, determine pathogenesis
        clinicalModel->update( *this, ageYears1, age0 == sim::zero() );
        clinicalModel->updateInfantDeaths( age0 );
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
    
    mon::reportMHI( mon::MHR_HOSTS, *this, 1 );
    mon::reportMHF( mon::MHF_AGE, *this, age(sim::now()).inYears() );
    bool patent = withinHostModel->summarize (*this);
    infIncidence->summarize (*this);
    
    if( patent ){
        // this should happen after all other reporting!
        removeFirstEvent( interventions::SubPopRemove::ON_FIRST_INFECTION );
    }
}

void Human::reportDeployment( ComponentId id, SimTime duration ){
    if( duration <= sim::zero() ) return; // nothing to do
    m_subPopExp[id] = sim::nowOrTs1() + duration;
    m_cohortSet = Survey::updateCohortSet( m_cohortSet, id, true );
}
void Human::removeFirstEvent( interventions::SubPopRemove::RemoveAtCode code ){
    const vector<ComponentId>& removeAtList = interventions::removeAtIds[code];
    for( vector<ComponentId>::const_iterator it = removeAtList.begin(), end = removeAtList.end(); it != end; ++it ){
        SubPopT::iterator expIt = m_subPopExp.find( *it );
        if( expIt != m_subPopExp.end() ){
            if( expIt->second > sim::nowOrTs0() ){
                // removeFirstEvent() is used for onFirstBout, onFirstTreatment
                // and onFirstInfection cohort options. Health system memory must
                // be reset for this to work properly; in theory the memory should
                // be independent for each cohort, but this is a usable approximation.
                flushReports();     // reset HS memory
                
                // report removal due to first infection/bout/treatment
                mon::reportMHI( mon::MHR_SUB_POP_REM_FIRST_EVENT, *this, 1 );
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
