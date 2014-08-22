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

#include "Clinical/ImmediateOutcomes.h"
#include "Clinical/CaseManagementCommon.h"
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
bool useDiagnosticUC = false;

double Params5Day::probGetsTreatment[Regimen::NUM];
double Params5Day::probParasitesCleared[Regimen::NUM-1];
double Params5Day::probClearedUcOnly;
double Params5Day::cureRateSevere;
WithinHost::TreatmentId Params5Day::treatments[Regimen::NUM];
ReportMeasureI measures[] = {
    Report::MI_TREATMENTS_1,    // first line official
    Report::MI_TREATMENTS_2,    // second line official
    Report::MI_TREATMENTS_1     // first line self treat
};


// ———  static, utility functions  ———

double getHealthSystemACRByName (const scnXml::TreatmentDetails& td, const string& drug)
{
    double val;
    if (drug == "CQ")
        val = td.getCQ().present() ? td.getCQ().get().getValue() : 0;
    else if (drug == "SP")
        val = td.getSP().present() ? td.getSP().get().getValue() : 0;
    else if (drug == "AQ")
        val = td.getAQ().present() ? td.getAQ().get().getValue() : 0;
    else if (drug == "SPAQ")
        val = td.getSPAQ().present() ? td.getSPAQ().get().getValue() : 0;
    else if (drug == "ACT")
        val = td.getACT().present() ? td.getACT().get().getValue() : 0;
    else if (drug == "QN")
        val = td.getQN().present() ? td.getQN().get().getValue() : 0;
    else if (drug == "selfTreatment")
        val = td.getSelfTreatment().getValue();
    else {
        throw util::xml_scenario_error ("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");
    }
    
    if( !(val >= 0.0 && val <= 1.0) ){
        throw util::xml_scenario_error (string("healthSystem initialACR/compliance/nonCompliersEffective: ").append(drug).append(" value must be in range [0,1]"));
    }
    return val;
}

// note: "CQOptional" is the same type as for other drugs
const scnXml::TreatmentActions::CQOptional& getHealthSystemTreatmentOptional(
    const scnXml::TreatmentActions& treatActs, const string& drug )
{
    if( drug == "CQ" ) return treatActs.getCQ();
    if( drug == "SP" ) return treatActs.getSP();
    if( drug == "AQ" ) return treatActs.getAQ();
    if( drug == "SPAQ" ) return treatActs.getSPAQ();
    if( drug == "ACT" ) return treatActs.getACT();
    if( drug == "QN" ) return treatActs.getQN();
    throw util::xml_scenario_error ("healthSystem.drugRegimen: firstLine / secondLine / inpatient has bad value");
}
WithinHost::TreatmentId getHealthSystemTreatmentByName( const scnXml::TreatmentActions& treatActs, const string& drug )
{
    const scnXml::TreatmentActions::CQOptional& elt = getHealthSystemTreatmentOptional( treatActs, drug );
    if( !elt.present() )
        throw util::xml_scenario_error( string("healthSystem.treatmentActions: description required for ").append(drug) );
    return WithinHost::WHInterface::addTreatment( elt.get() );
}


// ———  static, init  ———

void Params5Day::initParameters () {
    //TODO: make this compatible
    if (util::ModelOptions::option (util::INCLUDES_PK_PD))
        throw util::xml_scenario_error ("OldCaseManagement is not compatible with INCLUDES_PK_PD");
}

void Params5Day::setHealthSystem (const scnXml::HealthSystem& healthSystem) {
    if ( !healthSystem.getImmediateOutcomes().present() )
        throw util::xml_scenario_error ("Expected ImmediateOutcomes section in healthSystem data (initial or intervention)");
    const scnXml::HSImmediateOutcomes& hsioData = healthSystem.getImmediateOutcomes().get();
    
    const string &firstLine = hsioData.getDrugRegimen().getFirstLine(),
        &secondLine = hsioData.getDrugRegimen().getSecondLine(),
        &inpatient = hsioData.getDrugRegimen().getInpatient();
    
    const double pSeekOfficialCareUncomplicated1 = hsioData.getPSeekOfficialCareUncomplicated1().getValue();
    const double pSelfTreatment = hsioData.getPSelfTreatUncomplicated().getValue();
    

    // --- calculate probGetsTreatment ---

    probGetsTreatment[Regimen::SELF] = hsioData.getPSelfTreatUncomplicated().getValue();
    probGetsTreatment[Regimen::UC] = hsioData.getPSeekOfficialCareUncomplicated1().getValue() +
            hsioData.getPSelfTreatUncomplicated().getValue();
    probGetsTreatment[Regimen::UC2] = hsioData.getPSeekOfficialCareUncomplicated2().getValue();
    probGetsTreatment[Regimen::SEVERE] = hsioData.getPSeekOfficialCareSevere().getValue();
    if( !(
        pSeekOfficialCareUncomplicated1 >= 0.0 && pSelfTreatment >= 0.0
        && probGetsTreatment[Regimen::UC] <= 1.0
        && probGetsTreatment[Regimen::UC2] >= 0.0 && probGetsTreatment[Regimen::UC2] <= 1.0
        && probGetsTreatment[Regimen::SEVERE] >= 0.0 && probGetsTreatment[Regimen::SEVERE] <= 1.0
    ) ){
        throw util::xml_scenario_error ("healthSystem: pSeekOfficialCareXXX and pSelfTreatUncomplicated must be in range [0,1]");
    }

    // --- calculate probParasitesCleared ---

    const double complianceFirstLine = getHealthSystemACRByName (hsioData.getCompliance(), firstLine);
    const double complianceSecondLine = getHealthSystemACRByName (hsioData.getCompliance(), secondLine);

    const double cureRateFirstLine = getHealthSystemACRByName (hsioData.getInitialACR(), firstLine);
    const double cureRateSecondLine = getHealthSystemACRByName (hsioData.getInitialACR(), secondLine);

    const double nonCompliersEffectiveFirstLine = getHealthSystemACRByName (hsioData.getNonCompliersEffective(), firstLine);
    const double nonCompliersEffectiveSecondLine = getHealthSystemACRByName (hsioData.getNonCompliersEffective(), secondLine);

    const double complianceSelfTreatment = hsioData.getCompliance().getSelfTreatment().getValue();
    const double cureRateSelfTreatment = hsioData.getInitialACR().getSelfTreatment().getValue();
    if( !(
        complianceSelfTreatment >= 0.0 && complianceSelfTreatment <= 1.0
        && cureRateSelfTreatment >= 0.0 && cureRateSelfTreatment <= 1.0
    ) ){
        throw util::xml_scenario_error ("healthSystem initialACR/compliance/nonCompliersEffective: pSelfTreatment must be in range [0,1]");
    }
    
    
    //calculate probParasitesCleared 0
    if ( (pSeekOfficialCareUncomplicated1 + pSelfTreatment) > 0) {
        probParasitesCleared[Regimen::UC] = (pSeekOfficialCareUncomplicated1
                                   * (complianceFirstLine * cureRateFirstLine
                                      + (1 - complianceFirstLine) * nonCompliersEffectiveFirstLine)
                                   + pSelfTreatment
                                   * (complianceSelfTreatment * cureRateSelfTreatment
                                      + (1 - complianceSelfTreatment) * nonCompliersEffectiveFirstLine))
                                  / (pSeekOfficialCareUncomplicated1 + pSelfTreatment);
    } else {
        probParasitesCleared[Regimen::UC] = 0;
    }
    
    //calculate probParasitesCleared for uncomplicated cases
    probClearedUcOnly = complianceFirstLine * cureRateFirstLine
            + (1 - complianceFirstLine) * nonCompliersEffectiveFirstLine;
    probParasitesCleared[Regimen::UC2] = complianceSecondLine * cureRateSecondLine
            + (1 - complianceSecondLine) * nonCompliersEffectiveSecondLine;
    probParasitesCleared[Regimen::SELF] = complianceSelfTreatment * cureRateSelfTreatment
            + (1 - complianceSelfTreatment) * nonCompliersEffectiveFirstLine;
    
    cureRateSevere = getHealthSystemACRByName (hsioData.getInitialACR(), inpatient);

    treatments[Regimen::UC] = getHealthSystemTreatmentByName(hsioData.getTreatmentActions(), firstLine);
    if( secondLine == firstLine )
        treatments[Regimen::UC2] = treatments[Regimen::UC];
    else
        treatments[Regimen::UC2] = getHealthSystemTreatmentByName(hsioData.getTreatmentActions(), secondLine);
    treatments[Regimen::SELF] = treatments[Regimen::UC];
    if( inpatient == firstLine )
        treatments[Regimen::SEVERE] = treatments[Regimen::UC];
    else if( inpatient == secondLine )
        treatments[Regimen::SEVERE] = treatments[Regimen::UC2];
    else
        treatments[Regimen::SEVERE] = getHealthSystemTreatmentByName(hsioData.getTreatmentActions(), inpatient);
    
    useDiagnosticUC = hsioData.getUseDiagnosticUC();
    
    if( hsioData.getPrimaquine().present() ){
        if( !ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) )
            throw util::xml_scenario_error( "health-system's primaquine element only supported by vivax" );
        WithinHost::WHVivax::setHSParameters( hsioData.getPrimaquine().get() );
    }
}


// ———  per-human, construction and destruction  ———

ImmediateOutcomes::ImmediateOutcomes (double tSF) :
        _tLastTreatment (TimeStep::never),
        _treatmentSeekingFactor (tSF)
{}
ImmediateOutcomes::~ImmediateOutcomes() {
}


// ———  per-human, update  ———

void ImmediateOutcomes::doClinicalUpdate (Human& human, double ageYears) {
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
    
    if( _tLastTreatment == TimeStep::simulation ){
        human.removeFirstEvent( interventions::SubPopRemove::ON_FIRST_TREATMENT );
    }
    if( pgState & Episode::SICK ){
        human.removeFirstEvent( interventions::SubPopRemove::ON_FIRST_BOUT );
    }
}

void ImmediateOutcomes::uncomplicatedEvent (
    Human& human,
    Episode::State pgState
) {
    latestReport.update (human, Episode::State( pgState ) );

    Regimen::Type regimen = (_tLastTreatment + healthSystemMemory > TimeStep::simulation)
                            ? Regimen::UC2 : Regimen::UC ;
    
    double x = random::uniform_01();
    if( x < Params5Day::probGetsTreatment[regimen] * _treatmentSeekingFactor ){
        // UC1: official care OR self treatment
        // UC2: official care only
        
        if( useDiagnosticUC ){
            Survey::current().addInt( Report::MI_TREAT_DIAGNOSTICS, human, 1 );
            if( !human.withinHostModel->diagnosticDefault() )
                return; // negative outcome: no treatment
        }
        
        _tLastTreatment = TimeStep::simulation;
        Survey::current().addInt( measures[regimen], human, 1 );
        
        double y = random::uniform_01();
        if( y < Params5Day::probParasitesCleared[regimen] ){
            // Could report Episode::RECOVERY to latestReport,
            // but we don't report out-of-hospital recoveries anyway.
            human.withinHostModel->treatment( human, Params5Day::treatments[regimen] );
            if( regimen == Regimen::UC ){
                Survey::current().addInt( Report::MI_TREAT_SUCCESS_1, human, 1 );
            }
        } else {
            // No change in parasitological status: treated outside of hospital
        }
        if( regimen == Regimen::UC ){
            double p = x < Params5Day::probGetsTreatment[Regimen::SELF] ?
                Params5Day::probParasitesCleared[Regimen::SELF] :
                Params5Day::probClearedUcOnly;
            if( y < p ){
                Survey::current().addInt( Report::MI_TREAT_SUCCESS_2, human, 1 );
            }
        }
        
        if( human.withinHostModel->optionalPqTreatment() )
            Survey::current().addInt( Report::MI_PQ_TREATMENTS, human, 1 );
    } else {
        // No change in parasitological status: non-treated
    }
}

void ImmediateOutcomes::severeMalaria (
    Human &human,
    Episode::State pgState,
    double ageYears,
    int& doomed
) {
    double p2, p3, p4, p5, p6, p7;
    // Probability of getting treatment (only part which is case managment):
    p2 = Params5Day::probGetsTreatment[Regimen::SEVERE] * _treatmentSeekingFactor;
    // Probability of getting cured after getting treatment:
    p3 = Params5Day::cureRateSevere;
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
        _tLastTreatment = TimeStep::simulation;
        Survey::current().addInt( Report::MI_TREATMENTS_3, human, 1 );

        Episode::State stateTreated = Episode::State (pgState | Episode::EVENT_IN_HOSPITAL);
        if (q[5] <= prandom) { // Parasites cleared (treated, in hospital)
            human.withinHostModel->treatment( human, Params5Day::treatments[Regimen::SEVERE] );
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

void ImmediateOutcomes::massDrugAdministration( Human& human,
        Monitoring::ReportMeasureI screeningReport,
        Monitoring::ReportMeasureI drugReport )
{
    assert(false);      // should never be called
}

void ImmediateOutcomes::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    _tLastTreatment & stream;
    _treatmentSeekingFactor & stream;
}
void ImmediateOutcomes::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    _tLastTreatment & stream;
    _treatmentSeekingFactor & stream;
}

}
}
