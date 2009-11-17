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

#include "Clinical/Episode.h"
#include "Simulation.h"
#include "Surveys.h"

int Episode::reportingPeriodMemory;


Episode::Episode (istream& in) : _ageGroup(0) {
  in >> _time;
  if (_time != TIMESTEP_NEVER) {
    in >> _surveyPeriod;
    _ageGroup.read (in);
    int temp;
    in >> temp;
    _state = Pathogenesis::State (temp);
  }
}
Episode::~Episode ()
{
  report ();
}
void Episode::write (ostream& out) {
  out << _time << endl;
  if (_time != TIMESTEP_NEVER) {
    out << _surveyPeriod << endl;
    out << _ageGroup.i() << endl;
    out << _state << endl;
  }
}


void Episode::update (int simulationTime, SurveyAgeGroup ageGroup, Pathogenesis::State newState)
{
  if (simulationTime > (_time + reportingPeriodMemory)) {
    report ();
    
    _time = simulationTime;
    _surveyPeriod = Surveys.getSurveyPeriod();
    _ageGroup = ageGroup;
    _state = newState;
  } else {
    _state = Pathogenesis::State (_state | newState);
  }
}

void Episode::report () {
  if (_time == TIMESTEP_NEVER)	// Nothing to report
    return;
  
  if (_state & Pathogenesis::MALARIA) {
    if (_state & Pathogenesis::COMPLICATED)
      Surveys.at(_surveyPeriod)
	.reportSevereEpisodes (_ageGroup, 1);
    else // UC or UC2
      Surveys.at(_surveyPeriod)
	.reportUncomplicatedEpisodes (_ageGroup, 1);
  } else if (_state & Pathogenesis::SICK) {
    Surveys.at(_surveyPeriod)
      .reportNonMalariaFevers (_ageGroup, 1);
  }	// also possibility of nothing, but not reported in this case
  
  if (_state & Pathogenesis::EVENT_IN_HOSPITAL) {
    if (_state & Pathogenesis::DIRECT_DEATH) {
      Surveys.at(_surveyPeriod)
	.reportDirectDeaths (_ageGroup, 1)
	.reportHospitalDeaths (_ageGroup, 1);
    }
    else if (_state & Pathogenesis::SEQUELAE) {
      Surveys.at(_surveyPeriod)
	.reportSequelae (_ageGroup, 1)
	.reportHospitalSequelae (_ageGroup, 1);
    }
    else if (_state & Pathogenesis::RECOVERY)
      Surveys.at(_surveyPeriod)
	.reportHospitalRecoveries (_ageGroup, 1);
  } else {
    if (_state & Pathogenesis::DIRECT_DEATH)
      Surveys.at(_surveyPeriod)
	.reportDirectDeaths (_ageGroup, 1);
    else if (_state & Pathogenesis::SEQUELAE)
      Surveys.at(_surveyPeriod)
	.reportSequelae (_ageGroup, 1);
    // Don't care about out-of-hospital recoveries
  }
}
