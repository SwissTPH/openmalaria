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

#include "Clinical/EventScheduler.h"
#include "Clinical/ImmediateOutcomes.h"
#include "NeonatalMortality.h"
#include "Simulation.h"
#include "inputData.h"

int ClinicalModel::reportingPeriodMemory;
vector<int> ClinicalModel::infantIntervalsAtRisk;
vector<int> ClinicalModel::infantDeaths;


// -----  static methods  -----

void ClinicalModel::init () {
  infantDeaths.resize(Global::intervalsPerYear);
  infantIntervalsAtRisk.resize(Global::intervalsPerYear);
  PathogenesisModel::init();
  if (Global::modelVersion & CLINICAL_EVENT_SCHEDULER)
    ClinicalEventScheduler::init();
  else
    ClinicalImmediateOutcomes::initParameters();
}

ClinicalModel* ClinicalModel::createClinicalModel (double cF, double tSF) {
  if (Global::modelVersion & CLINICAL_EVENT_SCHEDULER)
    return new ClinicalEventScheduler (cF, tSF);
  else
    return new ClinicalImmediateOutcomes (cF, tSF);
}
ClinicalModel* ClinicalModel::createClinicalModel (istream& in) {
  if (Global::modelVersion & CLINICAL_EVENT_SCHEDULER)
    return new ClinicalEventScheduler (in);
  else
    return new ClinicalImmediateOutcomes (in);
}


// -----  non-static construction, destruction and checkpointing  -----

ClinicalModel::ClinicalModel (double cF) :
    pathogenesisModel(PathogenesisModel::createPathogenesisModel(cF)),
    _doomed(0)
{}
ClinicalModel::~ClinicalModel () {
  delete pathogenesisModel;
  // latestReport is reported, if any, by destructor
}

ClinicalModel::ClinicalModel (istream& in) :
    pathogenesisModel(PathogenesisModel::createPathogenesisModel(in))
{
  in >> latestReport;
  in >> _doomed; 
}
void ClinicalModel::write (ostream& out) {
  pathogenesisModel->write (out);
  out << latestReport << endl;
  out << _doomed << endl;
}


// -----  other non-static methods  -----

bool ClinicalModel::isDead (int ageTimeSteps) {
  if (ageTimeSteps > Global::maxAgeIntervals)	// too old (reached age limit)
    _doomed = DOOMED_TOO_OLD;
  if (_doomed > 0)	// killed by some means
    return true;	// remove from population
  return false;
}

void ClinicalModel::update (WithinHostModel& withinHostModel, double ageYears, int ageTimeSteps) {
  if (_doomed < 0)	// Countdown to indirect mortality
    _doomed -= Global::interval;
  
  //indirect death: if this human's about to die, don't worry about further episodes:
  if (_doomed <= -35) {	//clinical episode 6 intervals before
    latestReport.update(Simulation::simulationTime, Simulation::gMainSummary->ageGroup(ageYears), Pathogenesis::INDIRECT_MORTALITY, Outcome::INDIRECT_DEATH);
    _doomed = DOOMED_INDIRECT;
    return;
  }
  if(ageTimeSteps == 1) {
    // Chance of neonatal mortality:
    if (NeonatalMortality::eventNeonatalMortality()) {
      latestReport.update(Simulation::simulationTime, Simulation::gMainSummary->ageGroup(ageYears), Pathogenesis::INDIRECT_MORTALITY, Outcome::INDIRECT_DEATH);
      _doomed = DOOMED_NEONATAL;
      return;
    }
  }
  
  doClinicalUpdate (withinHostModel, ageYears);
}

void ClinicalModel::updateInfantDeaths (int ageTimeSteps) {
  // update array for the infant death rates
  if (ageTimeSteps <= (int)Global::intervalsPerYear){
    ++infantIntervalsAtRisk[ageTimeSteps-1];
    // Testing _doomed == -30 gives very slightly different results than
    // testing _doomed == DOOMED_INDIRECT (due to sim. end?)
    if (_doomed == DOOMED_COMPLICATED || _doomed == -30 || _doomed == DOOMED_NEONATAL){
      ++infantDeaths[ageTimeSteps-1];
    }
  }
}

void ClinicalModel::summarize (Summary& summary, double age) {
  pathogenesisModel->summarize (summary, age);
}
