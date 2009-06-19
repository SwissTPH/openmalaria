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

#include "Clinical/EventScheduler.h"
#include "inputData.h"
#include "GSLWrapper.h"
#include "WithinHostModel.h"
#include "simulation.h"


// -----  static init  -----

void ClinicalEventScheduler::init () {
  if (Global::interval != 1)
    throw xml_scenario_error("ClinicalEventScheduler is only designed for a 1-day timestep.");
  if (getCaseManagements() == NULL)
    throw xml_scenario_error ("ClinicalEventScheduler selected without caseManagements data in XML");
}


// -----  construction, destruction and checkpointing  -----

ClinicalEventScheduler::ClinicalEventScheduler (double cF, double tSF) :
    ClinicalModel (cF),
    pgState(Pathogenesis::NONE), reportState(Pathogenesis::NONE),
    pgChangeTimestep(TIMESTEP_NEVER), reportStartTimestep(TIMESTEP_NEVER)
{}
ClinicalEventScheduler::~ClinicalEventScheduler() {
}

ClinicalEventScheduler::ClinicalEventScheduler (istream& in) :
    ClinicalModel (in)
{
  int state;
  in >> state;
  pgState = (Pathogenesis::State) state;
  in >> state;
  reportState = (Pathogenesis::State) state;
  in >> pgChangeTimestep;
  in >> reportStartTimestep;
}
void ClinicalEventScheduler::write (ostream& out) {
  pathogenesisModel->write (out);
  out << latestReport;
  out << _doomed << endl; 
  out << pgState << endl;
  out << reportState << endl;
  out << pgChangeTimestep << endl;
  out << reportStartTimestep << endl;
}


// -----  other methods  -----

void ClinicalEventScheduler::doClinicalUpdate (WithinHostModel& withinHostModel, double ageYears) {
  // Run pathogenesisModel
  Pathogenesis::State newState = pathogenesisModel->determineState (ageYears, withinHostModel);
  
  /* Literally, if differences between (the combination of pgState and newState)
   * and pgState include SICK, MALARIA or COMPLICATED
   * or a second case of MALARIA and at least 3 days since last change */
  if ((((newState | pgState) ^ pgState) & (Pathogenesis::SICK | Pathogenesis::MALARIA | Pathogenesis::COMPLICATED)) ||
      (pgState & newState & Pathogenesis::MALARIA && Simulation::simulationTime >= pgChangeTimestep + 3))
  {
    if (Simulation::simulationTime > reportStartTimestep + 28) {	//NOTE: 28 is reporting length in days
      //FIXME: report
      reportState = pgState;
      reportStartTimestep = Simulation::simulationTime;
    }
    
    if ((newState & pgState) & Pathogenesis::MALARIA)
      newState = Pathogenesis::State (newState | Pathogenesis::SECOND_CASE);
    pgState = Pathogenesis::State (pgState | newState);
    pgChangeTimestep = Simulation::simulationTime;
    
    doCaseManagement (withinHostModel, ageYears);
  }
  
  
  if (pgState & Pathogenesis::INDIRECT_MORTALITY && _doomed == 0)
    _doomed = -Global::interval;	// start indirect mortality countdown
    
  
  if (pgState & Pathogenesis::COMPLICATED) {
    if (Simulation::simulationTime >= pgChangeTimestep + 10) {
      // FIXME: force recovery or something
      pgState = Pathogenesis::NONE;
    } else {
      //TODO: insert correct probabilities
      double rand = W_UNIFORM();
      if (rand < 0.5) {
	//TODO: if (rand < XXX && Simulation::simulationTime >= pgChangeTimestep + 5) report sequelae
	pgState = Pathogenesis::NONE;
      } else {
	if (rand < 0.75)
	  _doomed = 4;	// kill human (removed from simulation next timestep)
	// else stay in this state
      }
    }
  }
  
  // TODO: force leaving sick states after timeout
}

void ClinicalEventScheduler::doCaseManagement (WithinHostModel& withinHostModel, double ageYears) {
  //TODO: implement age-specificity of drug dosing
  //TODO: This is a rough and quick implementation, which could perhaps be improved.
  
  // Often individuals are not sick:
  if (!(pgState & Pathogenesis::SICK)) return;
  
  const scnXml::CaseManagements::CaseManagementSequence managements = getCaseManagements()->getCaseManagement();
  const scnXml::CaseManagement* caseManagement = NULL;
  for (scnXml::CaseManagements::CaseManagementConstIterator it = managements.begin(); it != managements.end(); ++it) {
    if (ageYears < it->getMaxAgeYrs() &&
      (!it->getMinAgeYrs().present() || it->getMinAgeYrs().get() <= ageYears))
      caseManagement = &*it;
  }
  if (caseManagement == NULL) {
    ostringstream msg;
    msg << "No case management for age " << ageYears;
    throw xml_scenario_error (msg.str());
  }
  
  const scnXml::CaseType::EndPointSequence* caseTypeSeq;
  if (pgState & Pathogenesis::MALARIA) {
    if (pgState & Pathogenesis::COMPLICATED)
      // FIXME: SEVERE/COINFECTION differences?
      caseTypeSeq = &caseManagement->getSev().getEndPoint();
    else if (pgState & Pathogenesis::SECOND_CASE)
      caseTypeSeq = &caseManagement->getUc2().getEndPoint();
    else
      caseTypeSeq = &caseManagement->getUc1().getEndPoint();
  } else /*if (pgState & Pathogenesis::SICK) [true by above check]*/ {	// sick but not from malaria
    caseTypeSeq = &caseManagement->getNmf().getEndPoint();
  }
  
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
  
  //FIXME: build list of decisions by ID and use (don't read every time)
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
  //TODO sort out reporting
  //latestReport.update(Simulation::simulationTime, agegroup, entrypoint, Outcome::PARASITES_PKPD_DEPENDENT_RECOVERS_OUTPATIENTS);
}
