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

#include "NewCaseManagement.h"
#include "WithinHostModel.h"
#include "inputData.h"
#include "GSLWrapper.h"


// -----  init  -----

void NewCaseManagement::init () {
  caseManagementMemory = get_health_system_memory();
  cerr << "Warning: NewCaseManagement has no way of determining clinical outcomes" << endl;
  if (!(Global::modelVersion & INCLUDES_PK_PD)) {
    throw xml_scenario_error ("Error: NewCaseManagement relies on INCLUDES_PK_PD to medicate treatment.");
  }
}

NewCaseManagement::NewCaseManagement(double tSF) :
    CaseManagementModel (tSF)
{
}

NewCaseManagement::NewCaseManagement(istream& in) :
  CaseManagementModel (in)
{
}

NewCaseManagement::~NewCaseManagement(){
}


void NewCaseManagement::doCaseManagement (Pathogenesis::State pgState, WithinHostModel& withinHostModel, Event&, double ageYears, int&) {
  //TODO: implement age-specificity of drug dosing
  //TODO: This is a rough and quick implementation, which could perhaps be improved.
  
  // Often individuals are not sick:
  if (pgState == Pathogenesis::NONE) return;
  
  //NOTE: should we just return in these cases? Maybe data should be read in init.
  if (getCaseManagements() == NULL) return;
  const scnXml::CaseManagements::CaseManagementSequence managements = getCaseManagements()->getCaseManagement();
  if (managements.size() == 0) return;
  const scnXml::CaseManagement* caseManagement = NULL;
  for (scnXml::CaseManagements::CaseManagementConstIterator it = managements.begin(); it != managements.end(); ++it)
    if (ageYears < it->getMaxAgeYrs() &&
      (!it->getMinAgeYrs().present() || it->getMinAgeYrs().get() <= ageYears))
      caseManagement = &*it;
    if (caseManagement == NULL) {
      ostringstream msg;
      msg << "No case management for age " << ageYears;
      throw xml_scenario_error (msg.str());
    }
  
  const scnXml::CaseType::EndPointSequence* caseTypeSeq;
  // FIXME: UC1/UC2 endpoints? (pgState & Pathogenesis::INDIRECT_MORTALITY)?
  if (pgState & Pathogenesis::MALARIA) {
    if (pgState & Pathogenesis::COMPLICATED)
      // FIXME: severe / copgState differences?
      caseTypeSeq = &caseManagement->getSev().getEndPoint();
    else //if (pgState & Pathogenesis::UNCOMPLICATED)
      caseTypeSeq = &caseManagement->getUc1().getEndPoint();
  } else if (pgState & Pathogenesis::NON_MALARIA) {
    caseTypeSeq = &caseManagement->getNmf().getEndPoint();
  } else {
    // NOTE: shouldn't happen; could just be checked in debug mode.
    cerr << "Invalid pgState code: "<<pgState<<endl;
    return;
  }
  /* UC2 should be the case sometimes:
    caseTypeSeq = &caseManagement->getUc2().getEndPoint();
  */
  
  double randCum = W_UNIFORM();
  int decisionID = -1;
  for (scnXml::CaseType::EndPointConstIterator it = caseTypeSeq->begin(); it != caseTypeSeq->end(); ++it) {
    randCum -= it->getP();
    if (randCum < 0) {
      decisionID = it->getDecision();
      break;
    }
  }
  if (decisionID < 0) {
    throw xml_scenario_error ("Sum of probabilities of case management end-points for some severity type less than 1");
  }
  
  //FIXME: build list of decisions by ID and use
  const scnXml::Decisions::DecisionSequence& decisions = caseManagement->getDecisions().getDecision();
  if ((int)decisions.size() < decisionID) {
    cerr << "A decision for a case-management end-point doesn't exist (number "<<decisionID<<")!" << endl;
    return;
  }
  const scnXml::Decision::MedicateSequence& medicates = decisions[decisionID-1].getMedicate();
  if (medicates.size() == 0)
    return;
  
  for (size_t medicateID=0;medicateID<medicates.size(); medicateID++) {
    double qty=medicates[medicateID].getQty();
    int time=medicates[medicateID].getTime();
    string name=medicates[medicateID].getName();
    withinHostModel.medicate(name,qty,time);
  }
}
