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

#include "Clinical/CaseManagementCommon.h"
#include "Clinical/EventScheduler.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Host/NeonatalMortality.h"

#include "inputData.h"
#include "Monitoring/Surveys.h"
#include "util/ModelOptions.h"

namespace OM { namespace Clinical {

vector<int> ClinicalModel::infantIntervalsAtRisk;
vector<int> ClinicalModel::infantDeaths;
double ClinicalModel::_nonMalariaMortality;


// -----  static methods  -----

void ClinicalModel::init () {
  infantDeaths.resize(TimeStep::stepsPerYear);
  infantIntervalsAtRisk.resize(TimeStep::stepsPerYear);
  _nonMalariaMortality=InputData.getParameter(Params::NON_MALARIA_INFANT_MORTALITY);
  
  Pathogenesis::PathogenesisModel::init();
  Episode::init();
  if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER))
    ClinicalEventScheduler::init();
  else
    ClinicalImmediateOutcomes::initParameters();
  CaseManagementCommon::initCommon();
}
void ClinicalModel::cleanup () {
  CaseManagementCommon::cleanupCommon();
    if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER))
	ClinicalEventScheduler::cleanup();
    Pathogenesis::PathogenesisModel::cleanup ();
}

void ClinicalModel::staticCheckpoint (istream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
    CaseManagementCommon::staticCheckpoint(stream);
}
void ClinicalModel::staticCheckpoint (ostream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
    CaseManagementCommon::staticCheckpoint(stream);
}

ClinicalModel* ClinicalModel::createClinicalModel (double cF, double tSF) {
  if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER))
    return new ClinicalEventScheduler (cF, tSF);
  else
    return new ClinicalImmediateOutcomes (cF, tSF);
}


void ClinicalModel::initMainSimulation () {
    for (TimeStep i(0);i<TimeStep::intervalsPerYear; ++i) {
	Clinical::ClinicalModel::infantIntervalsAtRisk[i.asInt()]=0;
	Clinical::ClinicalModel::infantDeaths[i.asInt()]=0;
    }
}

double ClinicalModel::infantAllCauseMort(){
  double infantPropSurviving=1.0;	// use to calculate proportion surviving
  for (TimeStep i(0);i<TimeStep::intervalsPerYear; ++i) {
    // multiply by proportion of infants surviving at each interval
    infantPropSurviving *= double(ClinicalModel::infantIntervalsAtRisk[i.asInt()]-ClinicalModel::infantDeaths[i.asInt()])
      / double(ClinicalModel::infantIntervalsAtRisk[i.asInt()]);
  }
  // Child deaths due to malaria (per 1000), plus non-malaria child deaths. Deaths per 1000 births is the return unit.
  return (1.0 - infantPropSurviving) * 1000.0 + _nonMalariaMortality;
}


// -----  non-static construction, destruction and checkpointing  -----

ClinicalModel::ClinicalModel (double cF) :
    pathogenesisModel(Pathogenesis::PathogenesisModel::createPathogenesisModel(cF)),
    _doomed(0)
{}
ClinicalModel::~ClinicalModel () {
  // latestReport is reported, if any, by destructor
}


// -----  other non-static methods  -----

bool ClinicalModel::isDead (TimeStep ageTimeSteps) {
  if (ageTimeSteps > TimeStep::maxAgeIntervals)	// too old (reached age limit)
    _doomed = DOOMED_TOO_OLD;
  if (_doomed > 0)	// killed by some means
    return true;	// remove from population
  return false;
}

void ClinicalModel::update (Human& human, double ageYears, TimeStep ageTimeSteps) {
  if (_doomed < 0)	// Countdown to indirect mortality
    _doomed -= TimeStep::interval;
  
  //indirect death: if this human's about to die, don't worry about further episodes:
  if (_doomed <= -35) {	//clinical bout 6 intervals before
     Monitoring::Surveys.getSurvey(human.getInCohort()).reportIndirectDeaths (human.getMonitoringAgeGroup(), 1);
    _doomed = DOOMED_INDIRECT;
    return;
  }
  if(ageTimeSteps == TimeStep(1)) {
    // Chance of neonatal mortality:
    if (Host::NeonatalMortality::eventNeonatalMortality()) {
      Monitoring::Surveys.getSurvey(human.getInCohort()).reportIndirectDeaths (human.getMonitoringAgeGroup(), 1);
      _doomed = DOOMED_NEONATAL;
      return;
    }
  }
  
  doClinicalUpdate (human, ageYears);
}

void ClinicalModel::updateInfantDeaths (TimeStep ageTimeSteps) {
  // update array for the infant death rates
  if (ageTimeSteps <= TimeStep::intervalsPerYear){
    ++infantIntervalsAtRisk[ageTimeSteps.asInt()-1];
    // Testing _doomed == -30 gives very slightly different results than
    // testing _doomed == DOOMED_INDIRECT (due to above if(..))
    if (_doomed == DOOMED_COMPLICATED || _doomed == -30 || _doomed == DOOMED_NEONATAL){
      ++infantDeaths[ageTimeSteps.asInt()-1];
    }
  }
}

void ClinicalModel::summarize (Monitoring::Survey& survey, Monitoring::AgeGroup ageGroup) {
  pathogenesisModel->summarize (survey, ageGroup);
}


void ClinicalModel::checkpoint (istream& stream) {
    (*pathogenesisModel) & stream;
    latestReport & stream;
    _doomed & stream;
}
void ClinicalModel::checkpoint (ostream& stream) {
    (*pathogenesisModel) & stream;
    latestReport & stream;
    _doomed & stream;
}

} }