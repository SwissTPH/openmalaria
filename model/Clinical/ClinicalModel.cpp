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

#include "Clinical/ClinicalModel.h"

#include "Clinical/CaseManagementCommon.h"
#include "Clinical/EventScheduler.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Host/NeonatalMortality.h"

#include "Monitoring/Surveys.h"
#include "util/ModelOptions.h"
#include <schema/scenario.h>

namespace OM { namespace Clinical {
    using namespace Monitoring;

vector<int> ClinicalModel::infantIntervalsAtRisk;
vector<int> ClinicalModel::infantDeaths;
double ClinicalModel::_nonMalariaMortality;

bool opt_event_scheduler = false;
bool ClinicalModel::indirectMortBugfix;

// -----  static methods  -----

void ClinicalModel::init( const Parameters& parameters, const scnXml::Model& model, const scnXml::HealthSystem& healthSystem ) {
    infantDeaths.resize(TimeStep::stepsPerYear);
    infantIntervalsAtRisk.resize(TimeStep::stepsPerYear);
    _nonMalariaMortality=parameters[Parameters::NON_MALARIA_INFANT_MORTALITY];
    indirectMortBugfix = util::ModelOptions::option (util::INDIRECT_MORTALITY_FIX);
    
    Episode::init( model.getClinical().getHealthSystemMemory() );
    if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER)){
        opt_event_scheduler = true;
        ClinicalEventScheduler::init( parameters, model );
    }else{
        ClinicalImmediateOutcomes::initParameters();
    }
    CaseManagementCommon::initCommon( parameters, healthSystem );
}
void ClinicalModel::cleanup () {
    CaseManagementCommon::cleanupCommon();
    if (opt_event_scheduler)
        ClinicalEventScheduler::cleanup();
}

void ClinicalModel::staticCheckpoint (istream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
}
void ClinicalModel::staticCheckpoint (ostream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
}

ClinicalModel* ClinicalModel::createClinicalModel (double tSF) {
  if (opt_event_scheduler)
    return new ClinicalEventScheduler (tSF);
  else
    return new ClinicalImmediateOutcomes (tSF);
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

ClinicalModel::ClinicalModel () :
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
        Surveys.getSurvey(human).addInt(
            Report::MI_INDIRECT_DEATHS, human.getMonitoringAgeGroup(), 1 );
        _doomed = DOOMED_INDIRECT;
        return;
    }
    if(ageTimeSteps == TimeStep(1) /* i.e. first update since birth */) {
        // Chance of neonatal mortality:
        if (Host::NeonatalMortality::eventNeonatalMortality()) {
            Surveys.getSurvey(human).addInt(
                Report::MI_INDIRECT_DEATHS, human.getMonitoringAgeGroup(), 1 );
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


void ClinicalModel::checkpoint (istream& stream) {
    latestReport & stream;
    _doomed & stream;
}
void ClinicalModel::checkpoint (ostream& stream) {
    latestReport & stream;
    _doomed & stream;
}

} }