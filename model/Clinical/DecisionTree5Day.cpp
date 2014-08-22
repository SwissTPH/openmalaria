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

enum CaseType { FirstLine, SecondLine, NumCaseTypes };

// These parameters are set by setHealthSystem() and do not need checkpointing.
double accessUCAny[NumCaseTypes];
double accessUCSelfTreat[NumCaseTypes];
double accessSevere;
auto_ptr<CMDecisionTree> treeUCOfficial, treeUCSelfTreat;
double cureRateSevere;
WithinHost::TreatmentId treatmentSevere;

// ———  static, init  ———

void DecisionTree5Day::setHealthSystem(const scnXml::HSDT5Day& hsDescription){
    double accessUCOfficial1 = hsDescription.getPSeekOfficialCareUncomplicated1().getValue();
    double accessUCOfficial2 = hsDescription.getPSeekOfficialCareUncomplicated2().getValue();
    // Note: this asymmetry is historical, and probably matters little:
    accessUCSelfTreat[FirstLine] = hsDescription.getPSelfTreatUncomplicated().getValue();
    accessUCSelfTreat[SecondLine] = 0.0;
    accessUCAny[FirstLine] = accessUCOfficial1 + accessUCSelfTreat[FirstLine];
    accessUCAny[SecondLine] = accessUCOfficial1 + accessUCSelfTreat[SecondLine];
    accessSevere = hsDescription.getPSeekOfficialCareSevere().getValue();
    
    if( accessUCOfficial1 < 0.0 || accessUCSelfTreat[FirstLine] < 0.0 || accessUCAny[FirstLine] > 1.0 ||
        accessUCOfficial2 < 0.0 || accessUCSelfTreat[SecondLine] < 0.0 || accessUCAny[SecondLine] > 1.0 ||
        accessSevere < 0.0 || accessSevere > 1.0 )
    {
        throw util::xml_scenario_error ("healthSystem: "
            "pSeekOfficialCareUncomplicated1 and pSelfTreatUncomplicated must be"
            " at least 0 and their sum must be at most 1, and"
            "pSeekOfficialCareUncomplicated2 and pSeekOfficialCareSevere must "
            "be in range [0,1]");
    }
    
    treeUCOfficial.reset( CMDecisionTree::create( hsDescription.getTreeUCOfficial() ).release() );
    treeUCSelfTreat.reset( CMDecisionTree::create( hsDescription.getTreeUCSelfTreat() ).release() );
    
    cureRateSevere = hsDescription.getCureRateSevere().getValue();
    treatmentSevere = WHInterface::addTreatment( hsDescription.getTreatmentSevere() );
    
    if( hsDescription.getPrimaquine().present() ){
        if( !ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) )
            throw util::xml_scenario_error( "health-system's primaquine element only supported by vivax" );
        WithinHost::WHVivax::setHSParameters( hsDescription.getPrimaquine().get() );
    }
}


// ———  per-human, construction and destruction  ———

DecisionTree5Day::DecisionTree5Day (double tSF) :
        m_tLastTreatment (TimeStep::never),
        m_treatmentSeekingFactor (tSF)
{}


// ———  per-human, update  ———

void DecisionTree5Day::doClinicalUpdate (Human& human, double ageYears) {
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

void DecisionTree5Day::uncomplicatedEvent ( Human& human, Episode::State pgState ){
    latestReport.update (human, Episode::State( pgState ) );
    
    // If last treatment prescribed was in recent memory, consider second line.
    CaseType regimen = FirstLine;
    if (m_tLastTreatment + healthSystemMemory > TimeStep::simulation){
        pgState = Episode::State (pgState | Episode::SECOND_CASE);
        regimen = SecondLine;
    }
    
    double x = random::uniform_01();
    if( x < accessUCAny[regimen] * m_treatmentSeekingFactor ){
        CMHostData hostData( human, human.getAgeInYears(), pgState );
        
        // Run tree (which may deploy treatment)
        CMDTOut output = ( x < accessUCSelfTreat[regimen] * m_treatmentSeekingFactor ) ?
            treeUCSelfTreat->exec( hostData ) :
            treeUCOfficial->exec( hostData );
        
        if( output.treated ){   // if any treatment or intervention deployed
            m_tLastTreatment = TimeStep::simulation;
            Monitoring::ReportMeasureI measure = (pgState & Episode::SECOND_CASE) ?
                Report::MI_TREATMENTS_1 : Report::MI_TREATMENTS_2;
            Survey::current().addInt( measure, human, 1 );
        }
        
        if( human.withinHostModel->optionalPqTreatment() )
            Survey::current().addInt( Report::MI_PQ_TREATMENTS, human, 1 );
    } else {
        // No care sought
    }
}

void DecisionTree5Day::severeMalaria (
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

void DecisionTree5Day::massDrugAdministration( Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport )
{
    assert(false);      // should never be called
}

void DecisionTree5Day::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    m_tLastTreatment & stream;
    m_treatmentSeekingFactor & stream;
}
void DecisionTree5Day::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    m_tLastTreatment & stream;
    m_treatmentSeekingFactor & stream;
}

}
}
