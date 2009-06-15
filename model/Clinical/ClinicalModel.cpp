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

#include "Clinical/ClinicalModel.h"

//#include "Clinical/EventScheduler.h"
#include "Clinical/ImmediateOutcomes.h"
#include "simulation.h"

// -----  static methods  -----

ClinicalModel* ClinicalModel::createClinicalModel (double cF, double tSF) {
  if (false)	//FIXME: condition
    return NULL;//new ClinicalEventScheduler (cF, tSF);
  else
    return new ClinicalImmediateOutcomes (cF, tSF);
}
ClinicalModel* ClinicalModel::createClinicalModel (istream& in) {
  if (false)
    return NULL;//new ClinicalEventScheduler (in);
  else
    return new ClinicalImmediateOutcomes (in);
}


// -----  non-static construction, destruction and checkpointing  -----

ClinicalModel::ClinicalModel (double cF, double tSF) :
    pathogenesisModel(PathogenesisModel::createPathogenesisModel(cF)),
    caseManagement(CaseManagementModel::createCaseManagementModel(tSF)),
    _doomed(0)
{}
ClinicalModel::~ClinicalModel () {
  delete pathogenesisModel;
  delete caseManagement; 
}

ClinicalModel::ClinicalModel (istream& in) :
    pathogenesisModel(PathogenesisModel::createPathogenesisModel(in)),
    caseManagement(CaseManagementModel::createCaseManagementModel(in))
{
  in >> latestReport;
  in >> _doomed; 
}
void ClinicalModel::write (ostream& out) {
  pathogenesisModel->write (out);
  caseManagement->write (out);
  out << latestReport;
  out << _doomed << endl; 
}


// -----  other non-static methods  -----

bool ClinicalModel::isDead (int ageTimeSteps) {
  if (ageTimeSteps > Global::maxAgeIntervals)	// too old (reached age limit)
    _doomed = 1;
  if (_doomed > 0)	// killed by some means
    return true;	// remove from population
  return false;
}

void ClinicalModel::update (WithinHostModel& withinHostModel, double ageYears, int ageTimeSteps) {
  if (_doomed < 0)	// Countdown to indirect mortality
    _doomed--;
  
  //indirect death: if this human's about to die, don't worry about further episodes:
  if (_doomed ==  -7) {	//clinical episode 6 intervals before
    latestReport.update(Simulation::simulationTime, Simulation::gMainSummary->ageGroup(ageYears), Diagnosis::INDIRECT_MALARIA_DEATH, Outcome::INDIRECT_DEATH);
    /*
    doomed=7 is the code for indirect death, and 6 for neonatal death.
    Individuals with positive doomed values are removed at the start of
    the next time step. They cannot be removed immediately because 
    their deaths need to be counted.
    */
    _doomed  = 7;
    //Indirect Neonatal mortality
    return;
  }
  // Neonatal mortality:
  if(ageTimeSteps == 1) {
    if (PathogenesisModel::eventNeonatalMortality()) {
      latestReport.update(Simulation::simulationTime, Simulation::gMainSummary->ageGroup(ageYears), Diagnosis::INDIRECT_MALARIA_DEATH, Outcome::INDIRECT_DEATH);
      _doomed  = 6;
      return;
    }
  }
  
  /* infectionEvent determines whether there is an acute episode, or concomitant fever and
  then whether the episode is severe, uncomplicated or there is an indirect
  death.
  doCaseManagement clears infections if there was an effective treatment, or calls medicate,
  and decides whether the patient lives, has sequelae, or dies.
  */
  caseManagement->doCaseManagement (pathogenesisModel->determineState (ageYears, withinHostModel),
				    withinHostModel,
				    latestReport,
				    ageYears,
				    _doomed);
}

void ClinicalModel::updateInfantDeaths (int ageTimeSteps) {
  // update array for the infant death rates
  if (ageTimeSteps <= (int)Global::intervalsPerYear){
    ++Global::infantIntervalsAtRisk[ageTimeSteps-1];
    if ((_doomed == 4) || (_doomed == -6) || (_doomed == 6)){
      ++Global::infantDeaths[ageTimeSteps-1];
    }
  }
}

void ClinicalModel::summarize (Summary& summary, double age) {
  pathogenesisModel->summarize (summary, age);
}
