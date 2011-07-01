/*
 This file is part of OpenMalaria.

 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "Clinical/ImmediateOutcomes.h"
#include "util/errors.h"
#include "util/ModelOptions.h"

namespace OM {
namespace Clinical {
using namespace ::OM::util;
using ::OM::Pathogenesis::State;

Diagnostic ClinicalImmediateOutcomes::massTreatDiagnostic;
double ClinicalImmediateOutcomes::probGetsTreatment[3];
double ClinicalImmediateOutcomes::probParasitesCleared[3];
double ClinicalImmediateOutcomes::cureRate[3];

// -----  static init  -----

void ClinicalImmediateOutcomes::initParameters () {
    if (util::ModelOptions::option (util::INCLUDES_PK_PD))
        throw util::xml_scenario_error ("OldCaseManagement is not compatible with INCLUDES_PK_PD");
}

void ClinicalImmediateOutcomes::setHealthSystem (const scnXml::HealthSystem& healthSystem) {
    if ( !healthSystem.getImmediateOutcomes().present() )
        throw util::xml_scenario_error ("Expected ImmediateOutcomes section in healthSystem data (initial or intervention)");

    setParasiteCaseParameters (healthSystem.getImmediateOutcomes().get());
}


// -----  construction and destruction  -----

ClinicalImmediateOutcomes::ClinicalImmediateOutcomes (double cF, double tSF) :
        ClinicalModel (cF),
        _tLastTreatment (TimeStep::never),
        _treatmentSeekingFactor (tSF)
{}
ClinicalImmediateOutcomes::~ClinicalImmediateOutcomes() {
}


// -----  other methods  -----

void ClinicalImmediateOutcomes::massDrugAdministration(Human& human) {
    if( !massTreatDiagnostic.isPositive( human.withinHostModel->getTotalDensity() ) ){
        return;
    }
    // We need to pass the is-severe state for the IPT code.
    human.withinHostModel->clearInfections(TimeStep::simulation, latestReport.getState() == Pathogenesis::STATE_SEVERE);
    Monitoring::Surveys.getSurvey(human.getInCohort()).reportMDA(human.getMonitoringAgeGroup(), 1);
}

void ClinicalImmediateOutcomes::doClinicalUpdate (Human& human, double ageYears) {
    bool effectiveTreatment = false;
    State pgState = pathogenesisModel->determineState (ageYears, *human.withinHostModel);

    if (pgState & Pathogenesis::MALARIA) {
        if (pgState & Pathogenesis::COMPLICATED)
            effectiveTreatment = severeMalaria (ageYears, human.getMonitoringAgeGroup(), _doomed, human.getInCohort());
        else if (pgState == Pathogenesis::STATE_MALARIA) {
            // NOTE: if condition means this doesn't happen if INDIRECT_MORTALITY is
            // included. Validity is debatable, but there's no point changing now.
            // (This does affect tests.)
            effectiveTreatment = uncomplicatedEvent (pgState, ageYears, human.getMonitoringAgeGroup(), human.getInCohort());
        }

        if ((pgState & Pathogenesis::INDIRECT_MORTALITY) && _doomed == 0)
            _doomed = -TimeStep::interval;

        if (util::ModelOptions::option (util::PENALISATION_EPISODES)) {
            human.withinHostModel->immunityPenalisation();
        }
    } else if (pgState & Pathogenesis::SICK) { // sick but not from malaria
        effectiveTreatment = uncomplicatedEvent (pgState, ageYears, human.getMonitoringAgeGroup(), human.getInCohort());
    }

    if (effectiveTreatment) {
        human.withinHostModel->clearInfections (TimeStep::simulation1(), latestReport.getState() == Pathogenesis::STATE_SEVERE);
    }

    if ( human.cohortFirstTreatmentOnly && _tLastTreatment == TimeStep::simulation ) {
        human.removeFromCohort();
    }
    if ( human.cohortFirstBoutOnly && (pgState & Pathogenesis::SICK) ) {
        human.removeFromCohort();
    }
}


// -----  private  -----

bool ClinicalImmediateOutcomes::uncomplicatedEvent (
    State pgState,
    double ageYears,
    Monitoring::AgeGroup ageGroup,
    bool inCohort
) {
    latestReport.update (inCohort, ageGroup,
                         State( pgState & Pathogenesis::STATE_MALARIA )   // mask to SICK and MALARIA flags
                        );

    Regimen::Type regimen = (_tLastTreatment + Episode::healthSystemMemory > TimeStep::simulation)
                            ? Regimen::UC2 : Regimen::UC
                            ;

    if ( probGetsTreatment[regimen]*_treatmentSeekingFactor > random::uniform_01() ) {
        _tLastTreatment = TimeStep::simulation;
        if ( regimen == Regimen::UC )
            Monitoring::Surveys.getSurvey(inCohort).reportTreatments1( ageGroup, 1 );
        if ( regimen == Regimen::UC2 )
            Monitoring::Surveys.getSurvey(inCohort).reportTreatments2( ageGroup, 1 );

        if (probParasitesCleared[regimen] > random::uniform_01()) {
            // Could report Pathogenesis::RECOVERY to latestReport,
            // but we don't report out-of-hospital recoveries anyway.
            return true;        // successful treatment
        } else {
            // No change in parasitological status: treated outside of hospital
            return false;
        }
    } else {
        // No change in parasitological status: non-treated
        return false;
    }
}

bool ClinicalImmediateOutcomes::severeMalaria (
    double ageYears,
    Monitoring::AgeGroup ageGroup,
    int& doomed,
    bool inCohort
) {
    Regimen::Type regimen = Regimen::SEVERE;

    double p2, p3, p4, p5, p6, p7;
    // Probability of getting treatment (only part which is case managment):
    p2 = probGetsTreatment[regimen] * _treatmentSeekingFactor;
    // Probability of getting cured after getting treatment:
    p3 = cureRate[regimen];
    // p4 is the hospital case-fatality rate from Tanzania
    p4 = caseFatality (ageYears);
    // p5 here is the community threshold case-fatality rate
    p5 = getCommunityCaseFatalityRate (p4);
    // p6 is P(seq) for treated patients
    p6 = pSequelaeInpatient (ageYears);
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
    q[3] = q[2] + p2 * p5 * (1 - p3);
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

    if (q[2] <= prandom) { // Patient gets in-hospital treatment
        _tLastTreatment = TimeStep::simulation;
        Monitoring::Surveys.getSurvey(inCohort).reportTreatments3( ageGroup, 1 );

        State sevTreated = State (Pathogenesis::STATE_SEVERE | Pathogenesis::EVENT_IN_HOSPITAL);
        if (q[5] <= prandom) { // Parasites cleared (treated, in hospital)
            if (q[6] > prandom) {
                latestReport.update (inCohort, ageGroup, State (sevTreated | Pathogenesis::DIRECT_DEATH));
                doomed  = 4;
            } else if (q[7] > prandom) { // Patient recovers, but with sequelae (don't report full recovery)
                latestReport.update (inCohort, ageGroup, State (sevTreated | Pathogenesis::SEQUELAE));
            } else { /*if (q[8] > prandom)*/
                latestReport.update (inCohort, ageGroup, State (sevTreated | Pathogenesis::RECOVERY));
            }
            return true;
        } else { // Treated but parasites not cleared (in hospital)
            if (q[3] > prandom) {
                latestReport.update (inCohort, ageGroup, State (sevTreated | Pathogenesis::DIRECT_DEATH));
                doomed  = 4;
            } else if (q[4] > prandom) { // sequelae without parasite clearance
                latestReport.update (inCohort, ageGroup, State (sevTreated | Pathogenesis::SEQUELAE));
            } else { /*if (q[5] > prandom)*/
                // No change in parasitological status: in-hospital patients
                latestReport.update (inCohort, ageGroup, Pathogenesis::STATE_SEVERE);
            }
            return false;
        }
    } else { // Not treated
        if (q[0] > prandom) {
            latestReport.update (inCohort, ageGroup, State (Pathogenesis::STATE_SEVERE | Pathogenesis::DIRECT_DEATH));
            doomed  = 4;
        } else if (q[1] > prandom) {
            latestReport.update (inCohort, ageGroup, State (Pathogenesis::STATE_SEVERE | Pathogenesis::SEQUELAE));
        } else { /*if (q[2] > prandom)*/
            // No change in parasitological status: non-treated
            latestReport.update (inCohort, ageGroup, Pathogenesis::STATE_SEVERE);
        }
        return false;
    }
}


double getHealthSystemACRByName (const scnXml::TreatmentDetails& td, string drug)
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
        throw util::xml_scenario_error ("healthSystem.drugRegimen->firstLine has bad value");
    }
    
    if( !(val >= 0.0 && val <= 1.0) ){
        throw util::xml_scenario_error (string("healthSystem initialACR/compliance/nonCompliersEffective: ").append(drug).append(" value must be in range [0,1]"));
    }
    return val;
}

void ClinicalImmediateOutcomes::setParasiteCaseParameters (const scnXml::HSImmediateOutcomes& hsioData)
{
    // --- calculate cureRate ---

    //We get the ACR depending on the name of firstLineDrug.
    cureRate[0] = getHealthSystemACRByName (hsioData.getInitialACR(),
                                            hsioData.getDrugRegimen().getFirstLine());

    //Calculate curerate 0
    const double pSeekOfficialCareUncomplicated1 = hsioData.getPSeekOfficialCareUncomplicated1().getValue();
    const double pSelfTreatment = hsioData.getPSelfTreatUncomplicated().getValue();
    if (pSeekOfficialCareUncomplicated1 + pSelfTreatment > 0) {
        double cureRateSelfTreatment = hsioData.getInitialACR().getSelfTreatment().getValue();

        cureRate[0] = (cureRate[0] * pSeekOfficialCareUncomplicated1
                       + cureRateSelfTreatment * pSelfTreatment)
                      / (pSeekOfficialCareUncomplicated1 + pSelfTreatment);
    }

    cureRate[1] = getHealthSystemACRByName (hsioData.getInitialACR(),
                                            hsioData.getDrugRegimen().getSecondLine());

    cureRate[2] = getHealthSystemACRByName (hsioData.getInitialACR(),
                                            hsioData.getDrugRegimen().getInpatient());


    // --- calculate probGetsTreatment ---

    probGetsTreatment[0] = hsioData.getPSeekOfficialCareUncomplicated1().getValue() + hsioData.getPSelfTreatUncomplicated().getValue();
    probGetsTreatment[1] = hsioData.getPSeekOfficialCareUncomplicated2().getValue();
    probGetsTreatment[2] = hsioData.getPSeekOfficialCareSevere().getValue();
    if( !(
        pSeekOfficialCareUncomplicated1 >= 0.0 && pSelfTreatment >= 0.0
        && probGetsTreatment[0] <= 1.0
        && probGetsTreatment[1] >= 0.0 && probGetsTreatment[1] <= 1.0
        && probGetsTreatment[2] >= 0.0 && probGetsTreatment[2] <= 1.0
    ) ){
        throw util::xml_scenario_error ("healthSystem: pSeekOfficialCareXXX and pSelfTreatUncomplicated must be in range [0,1]");
    }

    // --- calculate probParasitesCleared ---

    const string firstLineDrug = hsioData.getDrugRegimen().getFirstLine();
    const string secondLineDrug = hsioData.getDrugRegimen().getSecondLine();
    
    const double complianceFirstLine = getHealthSystemACRByName (hsioData.getCompliance(), firstLineDrug);
    const double complianceSecondLine = getHealthSystemACRByName (hsioData.getCompliance(), secondLineDrug);

    const double cureRateFirstLine = getHealthSystemACRByName (hsioData.getInitialACR(), firstLineDrug);
    const double cureRateSecondLine = getHealthSystemACRByName (hsioData.getInitialACR(), secondLineDrug);

    const double nonCompliersEffectiveFirstLine = getHealthSystemACRByName (hsioData.getNonCompliersEffective(), firstLineDrug);
    const double nonCompliersEffectiveSecondLine = getHealthSystemACRByName (hsioData.getNonCompliersEffective(), secondLineDrug);

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
        probParasitesCleared[0] = (pSeekOfficialCareUncomplicated1
                                   * (complianceFirstLine * cureRateFirstLine
                                      + (1 - complianceFirstLine) * nonCompliersEffectiveFirstLine)
                                   + pSelfTreatment
                                   * (complianceSelfTreatment * cureRateSelfTreatment
                                      + (1 - complianceSelfTreatment) * nonCompliersEffectiveFirstLine))
                                  / (pSeekOfficialCareUncomplicated1 + pSelfTreatment);
    } else {
        probParasitesCleared[0] = 0;
    }

    //calculate probParasitesCleared 1
    probParasitesCleared[1] = complianceSecondLine * cureRateSecondLine
                              + (1 - complianceSecondLine)
                              * nonCompliersEffectiveSecondLine;

    //calculate probParasitesCleared 2 : cool :)
    probParasitesCleared[2] = 0;
}


// -----  protected  -----

void ClinicalImmediateOutcomes::checkpoint (istream& stream) {
    ClinicalModel::checkpoint (stream);
    _tLastTreatment & stream;
    _treatmentSeekingFactor & stream;
}
void ClinicalImmediateOutcomes::checkpoint (ostream& stream) {
    ClinicalModel::checkpoint (stream);
    _tLastTreatment & stream;
    _treatmentSeekingFactor & stream;
}

}
}
