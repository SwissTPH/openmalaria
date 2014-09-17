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

#include "Clinical/DecisionTree5Day.h"
#include "Clinical/CaseManagementCommon.h"
#include "Clinical/CMDecisionTree.h"
#include "WithinHost/WHInterface.h"
#include "WithinHost/WHVivax.h"
#include "Monitoring/Survey.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/random.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;
using namespace Monitoring;

// These parameters are set by setHealthSystem() and do not need checkpointing.
Monitoring::ReportMeasureI CM5DayCommon::measures[NumCaseTypes] = {
    Report::MI_TREATMENTS_1,    // first line official
    Report::MI_TREATMENTS_2    // second line official
};
double CM5DayCommon::accessUCAny[NumCaseTypes];
double CM5DayCommon::accessUCSelfTreat[NumCaseTypes];
double CM5DayCommon::accessSevere;
double CM5DayCommon::cureRateSevere;
WithinHost::TreatmentId CM5DayCommon::treatmentSevere;


// ———  per-human, construction and destruction  ———

CM5DayCommon::CM5DayCommon (double tSF) :
        m_tLastTreatment (sim::never().ts()),
        m_treatmentSeekingFactor (tSF)
{}


// ———  per-human, update  ———

void CM5DayCommon::doClinicalUpdate (Human& human, double ageYears) {
    WithinHost::Pathogenesis::StatePair pg = human.withinHostModel->determineMorbidity( ageYears );
    Episode::State pgState = static_cast<Episode::State>( pg.state );

    if (pgState & Episode::MALARIA) {
        if (pgState & Episode::COMPLICATED){
            severeMalaria (human, pgState, ageYears, doomed);
        }else if (indirectMortBugfix || !pg.indirectMortality) {
            // NOTE: the "not indirect mortality" bit is a historical accident.
            // Validity is debatable, but there's no point changing now.
            // (This does affect tests.)
            uncomplicatedEvent (human, pgState);
        }

    } else if (pgState & Episode::SICK) { // sick but not from malaria
        uncomplicatedEvent (human, pgState);
    }
    
    if (pg.indirectMortality && doomed == NOT_DOOMED)
        doomed = -TimeStep::interval;
    
    if( m_tLastTreatment == TimeStep::simulation ){
        human.removeFirstEvent( interventions::SubPopRemove::ON_FIRST_TREATMENT );
    }
    if( pgState & Episode::SICK ){
        human.removeFirstEvent( interventions::SubPopRemove::ON_FIRST_BOUT );
    }
}

void CM5DayCommon::severeMalaria (
    Human &human,
    Episode::State pgState,
    double ageYears,
    int& doomed
) {
    double p2, p3, p4, p5, p6, p7;
    // Probability of getting treatment (only part which is case managment):
    p2 = accessSevere * m_treatmentSeekingFactor;
    // Probability of getting cured after getting treatment:
    p3 = cureRateSevere;
    // p4 is the hospital case-fatality rate from Tanzania
    p4 = caseFatalityRate.eval (ageYears);
    // p5 here is the community threshold case-fatality rate
    p5 = getCommunityCFR (p4);
    // p6 is P(seq) for treated patients
    p6 = pSequelaeInpatient.eval (ageYears);
    // p7 is P(seq) when parasites aren't cleared
    p7 = p6;

    double q[9];
    // Community deaths
    q[0] = (1 - p2) * p5;
    // Community sequelae
    q[1] = q[0] + (1 - p2) * (1 - p5) * p7;
    // Community survival
    q[2] = q[1] + (1 - p2) * (1 - p5) * (1 - p7);
    // Parasitological failure deaths
    q[3] = q[2] + p2 * (1 - p3) * p5;
    // Parasitological failure sequelae
    q[4] = q[3] + p2 * (1 - p3) * (1 - p5) * p7;
    // Parasitological failure survivors
    q[5] = q[4] + p2 * (1 - p3) * (1 - p5) * (1 - p7);
    // Parasitological success deaths
    q[6] = q[5] + p2 * p3 * p4;
    // Parasitological success sequelae
    q[7] = q[6] + p2 * p3 * (1 - p4) * p6;
    // Parasitological success survival
    q[8] = q[7] + p2 * p3 * (1 - p4) * (1 - p6);
    /*
    if (q(5).lt.1) stop
    NOT TREATED
    */

    double prandom = random::uniform_01();

    //NOTE: no diagnostics or PQ here; for now we only have severe when the patient dies
    if (q[2] <= prandom) { // Patient gets in-hospital treatment
        m_tLastTreatment = TimeStep::simulation;
        Survey::current().addInt( Report::MI_TREATMENTS_3, human, 1 );

        Episode::State stateTreated = Episode::State (pgState | Episode::EVENT_IN_HOSPITAL);
        if (q[5] <= prandom) { // Parasites cleared (treated, in hospital)
            human.withinHostModel->treatment( human, treatmentSevere );
            if (q[6] > prandom) {
                latestReport.update (human, Episode::State (stateTreated | Episode::DIRECT_DEATH));
                doomed  = DOOMED_COMPLICATED;
            } else if (q[7] > prandom) { // Patient recovers, but with sequelae (don't report full recovery)
                latestReport.update (human, Episode::State (stateTreated | Episode::SEQUELAE));
            } else { /*if (q[8] > prandom)*/
                latestReport.update (human, Episode::State (stateTreated | Episode::RECOVERY));
            }
        } else { // Treated but parasites not cleared (in hospital)
            if (q[3] > prandom) {
                latestReport.update (human, Episode::State (stateTreated | Episode::DIRECT_DEATH));
                doomed  = DOOMED_COMPLICATED;
            } else if (q[4] > prandom) { // sequelae without parasite clearance
                latestReport.update (human, Episode::State (stateTreated | Episode::SEQUELAE));
            } else { /*if (q[5] > prandom)*/
                // No change in parasitological status: in-hospital patients
                latestReport.update (human, pgState);
            }
        }
    } else { // Not treated
        if (q[0] > prandom) {
            latestReport.update (human, Episode::State (pgState | Episode::DIRECT_DEATH));
            doomed  = DOOMED_COMPLICATED;
        } else if (q[1] > prandom) {
            latestReport.update (human, Episode::State (pgState | Episode::SEQUELAE));
        } else { /*if (q[2] > prandom)*/
            // No change in parasitological status: non-treated
            latestReport.update (human, pgState);
        }
    }
}


// ———  per-human, intervention & checkpointing  ———

void CM5DayCommon::massDrugAdministration( Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport )
{
    assert(false);      // should never be called
}

void CM5DayCommon::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    m_tLastTreatment & stream;
    m_treatmentSeekingFactor & stream;
}
void CM5DayCommon::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    m_tLastTreatment & stream;
    m_treatmentSeekingFactor & stream;
}

}
}
