/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
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
#include "Host/HumanHet.h"
#include "Host/InfectionIncidenceModel.h"
#include "Clinical/ClinicalModel.h"
#include "Host/WithinHost/WHInterface.h"

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

// -----  Non-static functions: creation/destruction, checkpointing  -----

// Create new human
Human::Human(SimTime dateOfBirth) :
    infIncidence(InfectionIncidenceModel::createModel()),
    rng(util::master_RNG),
    dateOfBirth(dateOfBirth),
    toRemove(false),
    cohortSet(0),
    nextCtsDist(0)
{
    // Initial humans are created at time 0 and may have dateOfBirth in past. Otherwise dateOfBirth must be now.
    assert( dateOfBirth == sim::nowOrTs1() || (sim::now() == sim::zero() && dateOfBirth < sim::now()) );
    
    HumanHet het = HumanHet::sample(rng);
    withinHostModel = WithinHost::WHInterface::createWithinHostModel( rng, het.comorbidityFactor );
    auto iiFactor = infIncidence->getAvailabilityFactor(rng, 1.0);
    perHostTransmission.initialise (rng, het.availabilityFactor * iiFactor);
    clinicalModel = Clinical::ClinicalModel::createClinicalModel (het.treatmentSeekingFactor);
}

// -----  Non-static functions: per-time-step update  -----
vector<double> EIR_per_genotype;        // cache (not thread safe)

namespace human
{
    void checkpoint(Human &human, istream &stream)
    {
        human.perHostTransmission & stream;
        human.infIncidence & stream;
        human.withinHostModel & stream;
        human.clinicalModel & stream;
        human.rng.checkpoint(stream);
        human.dateOfBirth & stream;
        human.vaccine & stream;
        human.monitoringAgeGroup & stream;
        human.cohortSet & stream;
        human.nextCtsDist & stream;
        human.subPopExp & stream;
    }

    void checkpoint(Human &human, ostream &stream)
    {
        human.perHostTransmission & stream;
        human.infIncidence & stream;
        human.withinHostModel & stream;
        human.clinicalModel & stream;
        human.rng.checkpoint(stream);
        human.dateOfBirth & stream;
        human.vaccine & stream;
        human.monitoringAgeGroup & stream;
        human.cohortSet & stream;
        human.nextCtsDist & stream;
        human.subPopExp & stream;
    }

    void reportDeployment(Human &human, ComponentId id, SimTime duration )
    {
        if( duration <= sim::zero() ) return; // nothing to do
        human.subPopExp[id] = sim::nowOrTs1() + duration;
        human.cohortSet = mon::updateCohortSet( human.cohortSet, id, true );
    }

    void removeFirstEvent(Human &human, interventions::SubPopRemove::RemoveAtCode code )
    {
        const vector<ComponentId>& removeAtList = interventions::removeAtIds[code];
        for( auto it = removeAtList.begin(), end = removeAtList.end(); it != end; ++it ){
            auto expIt = human.subPopExp.find( *it );
            if( expIt != human.subPopExp.end() ){
                if( expIt->second > sim::nowOrTs0() ){
                    // removeFirstEvent() is used for onFirstBout, onFirstTreatment
                    // and onFirstInfection cohort options. Health system memory must
                    // be reset for this to work properly; in theory the memory should
                    // be independent for each cohort, but this is a usable approximation.
                    human.clinicalModel->flushReports();     // reset HS memory
                    
                    // report removal due to first infection/bout/treatment
                    mon::reportEventMHI( mon::MHR_SUB_POP_REM_FIRST_EVENT, human, 1 );
                }
                human.cohortSet = mon::updateCohortSet( human.cohortSet, expIt->first, false );
                // remove (affects reporting, restrictToSubPop and cumulative deployment):
                human.subPopExp.erase( expIt );
            }
        }
    }

    void summarize(Human &human, bool surveyOnlyNewEp) {
        if( surveyOnlyNewEp && human.clinicalModel->isExistingCase() ){
            // This modifies the denominator to treat the health-system-memory
            // period immediately after a bout as 'not at risk'.
            return;
        }
        
        mon::reportStatMHI( mon::MHR_HOSTS, human, 1 );
        mon::reportStatMHF( mon::MHF_AGE, human, sim::inYears(human.age(sim::now())) );
        bool patent = human.withinHostModel->summarize (human);
        human.infIncidence->summarize (human);
        
        if( patent && mon::isReported() ){
            // this should happen after all other reporting!
            human::removeFirstEvent(human, interventions::SubPopRemove::ON_FIRST_INFECTION);
        }
    }

    void update(Human &human, Transmission::TransmissionModel& transmission)
    {
        // For integer age checks we use age0 to e.g. get 73 steps comparing less than 1 year old
        SimTime age0 = human.age(sim::ts0());
        if (human.clinicalModel->isDead(age0)) {
            human.toRemove = true;
            return;
        }
        
        util::streamValidate( age0 );
        // Age at  the end of the update period. In most cases
        // the difference between this and age at the start is not especially
        // important in the model design, but since we parameterised with
        // ageYears1 we should stick with it.
        double ageYears1 = sim::inYears(human.age(sim::ts1()));
        // monitoringAgeGroup is the group for the start of the time step.
        human.monitoringAgeGroup.update( age0 );
        // check sub-pop expiry
        for( auto expIt = human.subPopExp.begin(), expEnd = human.subPopExp.end(); expIt != expEnd; ) {
            if( !(expIt->second >= sim::ts0()) ){       // membership expired
                // don't flush reports
                // report removal due to expiry
                mon::reportEventMHI( mon::MHR_SUB_POP_REM_TOO_OLD, human, 1 );
                human.cohortSet = mon::updateCohortSet( human.cohortSet, expIt->first, false );
                // erase element, but continue iteration
                expIt = human.subPopExp.erase( expIt );
            }else{
                ++expIt;
            }
        }
        // ageYears1 used only in PerHost::relativeAvailabilityAge(); difference to age0 should be minor
        double EIR = transmission.getEIR( human, age0, ageYears1, EIR_per_genotype );
        
        int nNewInfs = human.infIncidence->numNewInfections( human, EIR );
        
        // ageYears1 used when medicating drugs (small effect) and in immunity model (which was parameterised for it)
        human.withinHostModel->update(human, human.rng, nNewInfs, EIR_per_genotype, ageYears1);
        
        human.infIncidence->reportNumNewInfections(human, nNewInfs);
        
        // ageYears1 used to get case fatality and sequelae probabilities, determine pathogenesis
        human.clinicalModel->update( human, ageYears1, age0 == sim::zero() );
        human.clinicalModel->updateInfantDeaths( age0 );
    }
}

} }
