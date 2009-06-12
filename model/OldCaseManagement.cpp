/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include "OldCaseManagement.h"
#include "WithinHostModel.h"
#include "inputData.h"
#include "global.h"
#include "simulation.h"
#include "summary.h"
#include "GSLWrapper.h"
#include <limits>

const int OldCaseManagement::SEQUELAE_AGE_BOUND[NUM_SEQUELAE_AGE_GROUPS] = { 1, 10 };
double OldCaseManagement::_oddsRatioThreshold;
bool OldCaseManagement::_noMortality;
std::vector<double> OldCaseManagement::_inputAge;
std::vector<double> OldCaseManagement::_caseFatalityRate;
double OldCaseManagement::probGetsTreatment[3];
double OldCaseManagement::probParasitesCleared[3];
double OldCaseManagement::cureRate[3];
double OldCaseManagement::probSequelaeTreated[2];
double OldCaseManagement::probSequelaeUntreated[2];


// -----  init  -----

void OldCaseManagement::init (){
  if (Global::modelVersion & INCLUDES_PK_PD) {
    cerr << "Warning: OldCaseManagement's case management predetermines clinical outcomes, and is not currently compatible with INCLUDES_PK_PD" << endl;
  }
  _oddsRatioThreshold = exp(getParameter(Params::LOG_ODDS_RATIO_CF_COMMUNITY));
  
  setParasiteCaseParameters ();
  
  const scnXml::ByAgeItems::ItemSequence& items = getHealthSystem().getPSequelaeInpatient().getItem();
  for (int agegrp=0; agegrp < NUM_SEQUELAE_AGE_GROUPS; agegrp++) {
    for (size_t i = 0; i < items.size(); i++) {
      if (items[i].getMaxAgeYrs() > SEQUELAE_AGE_BOUND[agegrp]) {
        probSequelaeTreated[agegrp] =
        probSequelaeUntreated[agegrp] = items[i].getValue();
        goto gotItem;
      }
    }
    throw xml_scenario_error ("In scenario.xml: healthSystem: pSequelaeInpatient: expected item with maxAgeYrs > 10");
    gotItem:;	// alternative to for ... break ... else
  }
  
  readCaseFatalityRatio();
}

OldCaseManagement::OldCaseManagement(double tSF) :
    CaseManagementModel (tSF), _latestRegimen (0)
{
}

OldCaseManagement::OldCaseManagement(istream& in) :
    CaseManagementModel (in)
{
  in >> _latestRegimen; 
}

OldCaseManagement::~OldCaseManagement(){
}


void OldCaseManagement::doCaseManagement (Pathogenesis::Infection infection, WithinHostModel& withinHostModel, double ageYears, int& doomed) {
  bool effectiveTreatment =false;
  
  if (infection & Pathogenesis::MALARIA) {
    if (infection & Pathogenesis::COMPLICATED)
      effectiveTreatment=severeMalaria(ageYears, doomed);
    else if (infection == Pathogenesis::UNCOMPLICATED)
      effectiveTreatment=uncomplicatedEvent(true, ageYears);
    
    if ((infection & Pathogenesis::INDIRECT_MORTALITY) && doomed == 0)
      doomed=-1;
    
    if (Global::modelVersion & PENALISATION_EPISODES) {
      withinHostModel.immunityPenalisation();
    }
  } else if(infection & Pathogenesis::NON_MALARIA) {
    effectiveTreatment = uncomplicatedEvent(false, ageYears);
  }
  
  if (effectiveTreatment) {
    if (!(Global::modelVersion & INCLUDES_PK_PD))
      withinHostModel.clearInfections(_latestEvent.getDiagnosis() == Diagnosis::SEVERE_MALARIA);
  }
}

// -----  private  -----

bool OldCaseManagement::uncomplicatedEvent(bool isMalaria, double ageYears){
  //ageGroup is not optimized
  int agegroup=Simulation::gMainSummary->ageGroup(ageYears);
    int entrypoint = isMalaria ? Diagnosis::UNCOMPLICATED_MALARIA
                               : Diagnosis::NON_MALARIA_FEVER;
    int nextRegimen=getNextRegimen(Simulation::simulationTime, entrypoint, _tLastTreatment, _latestRegimen);
    if (probGetsTreatment[nextRegimen-1]*_treatmentSeekingFactor > (W_UNIFORM())){
      _latestRegimen=nextRegimen;
      _tLastTreatment=Simulation::simulationTime;
      Simulation::gMainSummary->reportTreatment(agegroup, _latestRegimen);
      
      if (Global::modelVersion & INCLUDES_PK_PD){
        _latestEvent.update(Simulation::simulationTime, agegroup, entrypoint, Outcome::PARASITES_PKPD_DEPENDENT_RECOVERS_OUTPATIENTS);
        /*
        TODO: uncomplicatedEvent forces a call of pk PD model
        in the event that there is no treatment it should remain .false.
          TODO: 
        call medicate(currInd, "cq ", 300.0, 0)
        */
        return true;
      }
      else {
        if (probParasitesCleared[nextRegimen-1] > W_UNIFORM()){
          _latestEvent.update(Simulation::simulationTime, agegroup, entrypoint, Outcome::PARASITES_ARE_CLEARED_PATIENT_RECOVERS_OUTPATIENTS);
          return true;
        }
        else {
          _latestEvent.update(Simulation::simulationTime, agegroup, entrypoint, Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_OUTPATIENTS);
        }
      }
    }
    else {
      _latestEvent.update(Simulation::simulationTime, agegroup, entrypoint, Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_NON_TREATED);
    }
  return false;
}

bool OldCaseManagement::severeMalaria(double ageYears, int& doomed){
  /*
  DOCU
  Set doomed=4 if the patient dies.
  */

  int agegroup=Simulation::gMainSummary->ageGroup(ageYears);
  int isAdultIndex=1;
  if (ageYears >= 5.0) {
    isAdultIndex=0;
  }
  /*
  TODO: is there a reason we cannot first decide if the patient is treated, then increase latestRegimen
  and latestTreatment, instead of resetting it if not treated? (Can't think of one, but
  do we want to change this section of code rather than just introducing the new alternative (TS))
  */ 
  int nextRegimen=getNextRegimen(Simulation::simulationTime, Diagnosis::SEVERE_MALARIA, _tLastTreatment, _latestRegimen);
  
  double p2, p3, p4, p5, p6, p7;
  // Probability of getting treatment (only part which is case managment):
  p2=probGetsTreatment[nextRegimen-1]*_treatmentSeekingFactor;
  // Probability of getting cured after getting treatment:
  p3=cureRate[nextRegimen-1];
  // p4 is the hospital case-fatality rate from Tanzania
  p4=caseFatality(ageYears);
  // p5 here is the community threshold case-fatality rate
  p5=getCommunityCaseFatalityRate(p4);
  p6=probSequelaeTreated[isAdultIndex];
  p7=probSequelaeUntreated[isAdultIndex];
  
  double q[9];
  //	Community deaths
  q[0]=(1-p2)*p5;
  //	Community sequelae
  q[1]=q[0]+(1-p2)*(1-p5)*p7;
  //	Community survival
  q[2]=q[1]+(1-p2)*(1-p5)*(1-p7);
  //	Parasitological failure deaths
  q[3]=q[2]+p2*p5*(1-p3);
  //	Parasitological failure sequelae
  q[4]=q[3]+p2*(1-p3)*(1-p5)*p7;
  //	Parasitological failure survivors
  q[5]=q[4]+p2*(1-p3)*(1-p5)*(1-p7);
  //	Parasitological success deaths
  q[6]=q[5]+p2*p3*p4;
  //	Parasitological success sequelae
  q[7]=q[6]+p2*p3*(1-p4)*p6;
  //	Parasitological success survival
  q[8]=q[7]+p2*p3*(1-p4)*(1-p6);
  /*
  if (q(5).lt.1) stop
  NOT TREATED
  */
  
  double prandom=(W_UNIFORM());
    
  if ( q[0] >  prandom) {
    _latestEvent.update(Simulation::simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PATIENT_DIES_NON_TREATED);
    doomed  = 4;
  }
  else if( q[1] >  prandom) {
    _latestEvent.update(Simulation::simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_NOT_CLEARED_PATIENT_HAS_SEQUELAE_NON_TREATED);
  }
  else if( q[2] >  prandom) {
    _latestEvent.update(Simulation::simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_NON_TREATED);
  }
  else {
    _tLastTreatment = Simulation::simulationTime;
    _latestRegimen = nextRegimen;
    Simulation::gMainSummary->reportTreatment(agegroup,_latestRegimen);
    
    if( q[3] >  prandom) {
      _latestEvent.update(Simulation::simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PATIENT_DIES_INPATIENTS);
      doomed  = 4;
    }
    else if( q[4] >  prandom) {
      _latestEvent.update(Simulation::simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_NOT_CLEARED_PATIENT_HAS_SEQUELAE_INPATIENTS);
    }
    else if( q[5] >  prandom) {
      _latestEvent.update(Simulation::simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_INPATIENTS);
    }
    else {
      if( q[6] >  prandom) {
        _latestEvent.update(Simulation::simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PATIENT_DIES_INPATIENTS);
        doomed  = 4;
      }
      else if( q[7] >  prandom) {
        _latestEvent.update(Simulation::simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_ARE_CLEARED_PATIENT_HAS_SEQUELAE_INPATIENTS);
      }
      else // assume true, so we don't get another else case (DH): if( q[8] >=  prandom)
      {
        _latestEvent.update(Simulation::simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_ARE_CLEARED_PATIENT_RECOVERS_INPATIENTS);
      }
      return true;
    }
  }
  return false;
}

void OldCaseManagement::readCaseFatalityRatio(){
  const scnXml::AgeGroups::GroupSequence& xmlGroupSeq = getHealthSystem().getCFR().getGroup();
  
  int numOfGroups = xmlGroupSeq.size();
  _inputAge.resize(numOfGroups+1);
  _caseFatalityRate.resize(numOfGroups+1);
  
  for (int i=0;i<numOfGroups; i++) {
    const scnXml::Group& xmlGroup = xmlGroupSeq[i];
    _inputAge[i] = xmlGroup.getLowerbound();
    _caseFatalityRate[i] = xmlGroup.getCfr();
  }
  
  //NOTE: need to make sure _inputAge[0] == 0 (or less)
  _inputAge[0] = 0;
  
  //Use a high number for the upper bound.
  _inputAge[numOfGroups] = numeric_limits<double>::infinity();
  //CFR is constant for everyone above the highest upperbound
  _caseFatalityRate[numOfGroups] = _caseFatalityRate[numOfGroups-1];
  _noMortality = (numOfGroups == 1) && ( _caseFatalityRate[0] == 0);
}

double OldCaseManagement::getCommunityCaseFatalityRate(double caseFatalityRatio) {
  double x = caseFatalityRatio * _oddsRatioThreshold;
  return x / (1 - caseFatalityRatio + x);
}

int OldCaseManagement::getNextRegimen(int simulationTime, int diagnosis, int tLastTreated, int regimen) {
  if (diagnosis == Diagnosis::SEVERE_MALARIA)
    return 3;
  
  if (tLastTreated > (simulationTime-caseManagementMemory))
    return 2;
  
  return 1;
}

double OldCaseManagement::caseFatality(double ageYears) {
  // NOTE: assume ageYears >= 0 and _inputAge[0] <= 0
  if (_noMortality)
    return 0.0;
  
  size_t i = 0;
  while(_inputAge[i] <= ageYears)
    ++i;
  
  // _inputAge[i-1] <= ageYears < _inputAge[i]
  double a0 = _inputAge[i-1];
  double f0 = _caseFatalityRate[i-1];
  return (ageYears - a0) / (_inputAge[i] - a0) * (_caseFatalityRate[i]-f0) + f0;
}


double getHealthSystemACRByName (const scnXml::TreatmentDetails& td, string firstLineDrug) {
  if (firstLineDrug == "CQ")
    return td.getCQ().present() ? td.getCQ().get().getValue() : 0;
  else if (firstLineDrug == "SP")
    return td.getSP().present() ? td.getSP().get().getValue() : 0;
  else if (firstLineDrug == "AQ")
    return td.getAQ().present() ? td.getAQ().get().getValue() : 0;
  else if (firstLineDrug == "SPAQ")
    return td.getSPAQ().present() ? td.getSPAQ().get().getValue() : 0;
  else if (firstLineDrug == "ACT")
    return td.getACT().present() ? td.getACT().get().getValue() : 0;
  else if (firstLineDrug == "QN")
    return td.getQN().present() ? td.getQN().get().getValue() : 0;
  else if (firstLineDrug == "selfTreatment")
    return td.getSelfTreatment().getValue();
  else {
    throw xml_scenario_error ("healthSystem.drugRegimen->firstLine has bad value");
  }
}

void OldCaseManagement::setParasiteCaseParameters () {
  const scnXml::HealthSystem& healthSystem = getHealthSystem();
  
  
  // --- calculate cureRate ---
  
  //We get the ACR depending on the name of firstLineDrug.
  cureRate[0] = getHealthSystemACRByName (healthSystem.getInitialACR(),
                                          healthSystem.getDrugRegimen().getFirstLine());
  
  //Calculate curerate 0
  double pSeekOfficialCareUncomplicated1 = healthSystem.getPSeekOfficialCareUncomplicated1().getValue();
  double pSelfTreatment = healthSystem.getPSelfTreatUncomplicated().getValue();
  if (pSeekOfficialCareUncomplicated1 + pSelfTreatment > 0){
    double cureRateSelfTreatment = healthSystem.getInitialACR().getSelfTreatment().getValue();
    
    cureRate[0] = (cureRate[0] * pSeekOfficialCareUncomplicated1
        + cureRateSelfTreatment * pSelfTreatment)
       / (pSeekOfficialCareUncomplicated1+pSelfTreatment);
  }

  cureRate[1] = getHealthSystemACRByName (healthSystem.getInitialACR(),
                                          healthSystem.getDrugRegimen().getSecondLine());

  cureRate[2] = getHealthSystemACRByName (healthSystem.getInitialACR(),
                                          healthSystem.getDrugRegimen().getInpatient());
  
  
  // --- calculate probGetsTreatment ---
  
  probGetsTreatment[0] = healthSystem.getPSeekOfficialCareUncomplicated1().getValue() + healthSystem.getPSelfTreatUncomplicated().getValue();
  probGetsTreatment[1] = healthSystem.getPSeekOfficialCareUncomplicated2().getValue();
  probGetsTreatment[2] = healthSystem.getPSeekOfficialCareSevere().getValue();
  
  
  // --- calculate probParasitesCleared ---
  
  string firstLineDrug = healthSystem.getDrugRegimen().getFirstLine();
  string secondLineDrug = healthSystem.getDrugRegimen().getSecondLine();	
  pSeekOfficialCareUncomplicated1 = healthSystem.getPSeekOfficialCareUncomplicated1().getValue();

  double complianceFirstLine = getHealthSystemACRByName (healthSystem.getCompliance(), firstLineDrug);
  double complianceSecondLine = getHealthSystemACRByName (healthSystem.getCompliance(), secondLineDrug);

  double cureRateFirstLine = getHealthSystemACRByName (healthSystem.getInitialACR(), firstLineDrug);
  double cureRateSecondLine = getHealthSystemACRByName (healthSystem.getInitialACR(), secondLineDrug);

  double nonCompliersEffectiveFirstLine = getHealthSystemACRByName (healthSystem.getNonCompliersEffective(), firstLineDrug);
  double nonCompliersEffectiveSecondLine = getHealthSystemACRByName (healthSystem.getNonCompliersEffective(), secondLineDrug);

  pSelfTreatment = healthSystem.getPSelfTreatUncomplicated().getValue();
  double complianceSelfTreatment = healthSystem.getCompliance().getSelfTreatment().getValue();
  double cureRateSelfTreatment = healthSystem.getInitialACR().getSelfTreatment().getValue();
  
  //calculate probParasitesCleared 0
  if ((pSeekOfficialCareUncomplicated1 + pSelfTreatment) > 0){
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


// -----  checkpointing  -----

void OldCaseManagement::write(ostream& out) const {
  out << _latestEvent;
  out << _treatmentSeekingFactor << endl; 
  out << _tLastTreatment << endl; 
  out << _latestRegimen << endl; 
}
