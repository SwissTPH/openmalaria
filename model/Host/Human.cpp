/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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
#include "util/ModelOptions.h"
#include "util/vectors.h"
#include "util/StreamValidator.h"
#include "Population.h"
#include "interventions/InterventionManager.hpp"
#include "mon/reporting.h"
#include "schema/scenario.h"

namespace OM { namespace Host {
    using namespace OM::util;
    using interventions::ComponentId;
    
    bool surveyOnlyNewEp = false;

// -----  Static functions  -----

void Human::init( const Parameters& parameters, const scnXml::Scenario& scenario ){    // static
    HumanHet::init();
    surveyOnlyNewEp = scenario.getMonitoring().getSurveyOptions().getOnlyNewEpisode();
    
    const scnXml::Model& model = scenario.getModel();
    // Init models used by humans:
    Transmission::PerHost::init( model.getHuman().getAvailabilityToMosquitoes() );
    InfectionIncidenceModel::init( parameters );
    WithinHost::WHInterface::init( parameters, scenario );
    Clinical::ClinicalModel::init( parameters, scenario );
}


// -----  Non-static functions: creation/destruction, checkpointing  -----

// Create new human
Human::Human(uint64_t seed1, uint64_t seed2, SimTime dateOfBirth) :
    infIncidence(InfectionIncidenceModel::createModel()),
    m_rng(seed1, seed2),
    m_DOB(dateOfBirth),
    m_cohortSet(0),
    nextCtsDist(0)
{
    // Initial humans are created at time 0 and may have DOB in past. Otherwise DOB must be now.
    assert( m_DOB == sim::nowOrTs1() || (sim::now() == SimTime::zero() && m_DOB < sim::now()) );
    
    HumanHet het = HumanHet::sample(m_rng);
    withinHostModel = WithinHost::WHInterface::createWithinHostModel( m_rng, het.comorbidityFactor );
    auto iiFactor = infIncidence->getAvailabilityFactor(m_rng, 1.0);
    perHostTransmission.initialise (m_rng, het.availabilityFactor * iiFactor);
    clinicalModel = Clinical::ClinicalModel::createClinicalModel (het.treatmentSeekingFactor);
}

Human::Human(SimTime dateOfBirth, int dummy) :
    withinHostModel(nullptr),
    infIncidence(nullptr),
    clinicalModel(nullptr),
    m_rng(0, 0),
    m_DOB(dateOfBirth),
    m_cohortSet(0),
    nextCtsDist(0)
{}


// -----  Non-static functions: per-time-step update  -----

vector<double> EIR_per_genotype;        // cache (not thread safe)

bool Human::update(const Transmission::TransmissionModel& transmission, bool doUpdate) {
    // For integer age checks we use age0 to e.g. get 73 steps comparing less than 1 year old
    SimTime age0 = age(sim::ts0());
    if (clinicalModel->isDead(age0))
        return true;
    
    if (doUpdate){
        util::streamValidate( age0.inDays() );
        // Age at  the end of the update period. In most cases
        // the difference between this and age at the start is not especially
        // important in the model design, but since we parameterised with
        // ageYears1 we should stick with it.
        double ageYears1 = age(sim::ts1()).inYears();
        // monitoringAgeGroup is the group for the start of the time step.
        monitoringAgeGroup.update( age0 );
        // check sub-pop expiry
        for( auto expIt = m_subPopExp.begin(), expEnd = m_subPopExp.end(); expIt != expEnd; ) {
            if( !(expIt->second >= sim::ts0()) ){       // membership expired
                // don't flush reports
                // report removal due to expiry
                mon::reportEventMHI( mon::MHR_SUB_POP_REM_TOO_OLD, *this, 1 );
                m_cohortSet = mon::updateCohortSet( m_cohortSet, expIt->first, false );
                // erase element, but continue iteration
                expIt = m_subPopExp.erase( expIt );
            }else{
                ++expIt;
            }
        }
        // ageYears1 used only in PerHost::relativeAvailabilityAge(); difference to age0 should be minor
        double EIR = transmission.getEIR( *this, age0, ageYears1,
                EIR_per_genotype );
        int nNewInfs = infIncidence->numNewInfections( *this, EIR );
        
        // ageYears1 used when medicating drugs (small effect) and in immunity model (which was parameterised for it)
        withinHostModel->update(m_rng, nNewInfs, EIR_per_genotype, ageYears1,
                _vaccine.getFactor(interventions::Vaccine::BSV));
        
        // ageYears1 used to get case fatality and sequelae probabilities, determine pathogenesis
        clinicalModel->update( *this, ageYears1, age0 == SimTime::zero() );
        clinicalModel->updateInfantDeaths( age0 );
    }
    return false;
}

void Human::addInfection(){
    withinHostModel->importInfection(m_rng);
}

void Human::clearImmunity(){
    withinHostModel->clearImmunity();
}


void Human::summarize() {
    if( surveyOnlyNewEp && clinicalModel->isExistingCase() ){
        // This modifies the denominator to treat the health-system-memory
        // period immediately after a bout as 'not at risk'.
        return;
    }
    
    mon::reportStatMHI( mon::MHR_HOSTS, *this, 1 );
    mon::reportStatMHF( mon::MHF_AGE, *this, age(sim::now()).inYears() );
    bool patent = withinHostModel->summarize (*this);
    infIncidence->summarize (*this);
    
    if( patent && mon::isReported() ){
        // this should happen after all other reporting!
        removeFirstEvent( interventions::SubPopRemove::ON_FIRST_INFECTION );
    }
}

void Human::reportDeployment( ComponentId id, SimTime duration ){
    if( duration <= SimTime::zero() ) return; // nothing to do
    m_subPopExp[id] = sim::nowOrTs1() + duration;
    m_cohortSet = mon::updateCohortSet( m_cohortSet, id, true );
}
void Human::removeFirstEvent( interventions::SubPopRemove::RemoveAtCode code ){
    const vector<ComponentId>& removeAtList = interventions::removeAtIds[code];
    for( auto it = removeAtList.begin(), end = removeAtList.end(); it != end; ++it ){
        auto expIt = m_subPopExp.find( *it );
        if( expIt != m_subPopExp.end() ){
            if( expIt->second > sim::nowOrTs0() ){
                // removeFirstEvent() is used for onFirstBout, onFirstTreatment
                // and onFirstInfection cohort options. Health system memory must
                // be reset for this to work properly; in theory the memory should
                // be independent for each cohort, but this is a usable approximation.
                flushReports();     // reset HS memory
                
                // report removal due to first infection/bout/treatment
                mon::reportEventMHI( mon::MHR_SUB_POP_REM_FIRST_EVENT, *this, 1 );
            }
            m_cohortSet = mon::updateCohortSet( m_cohortSet, expIt->first, false );
            // remove (affects reporting, restrictToSubPop and cumulative deployment):
            m_subPopExp.erase( expIt );
        }
    }
}

void Human::flushReports (){
    clinicalModel->flushReports();
}

} }
