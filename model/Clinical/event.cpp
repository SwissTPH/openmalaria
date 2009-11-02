/*

 This file is part of OpenMalaria.

 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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
#include "Clinical/event.h"
#include "Simulation.h"
#include "inputData.h"
#include "util/gsl.h"
#include <algorithm>
#include "Clinical/ClinicalModel.h"
#include "summary.h"

void Event::update(int simulationTime, int ageGroup, Pathogenesis::State newState, int outcome){
  if ((newState == Pathogenesis::INDIRECT_MORTALITY) || (simulationTime>(_time + ClinicalModel::reportingPeriodMemory))){
    if (_time!=TIMESTEP_NEVER){
      Simulation::gMainSummary->report(*this);
    }
    _time=simulationTime;
    _surveyPeriod=Simulation::gMainSummary->getSurveyPeriod();
    _ageGroup=ageGroup;
    _state=newState;
    _outcome=outcome;
  }
  else {
    _outcome=std::max(outcome, _outcome);
    _state = Pathogenesis::State(_state | newState);
  }
}

Event::~Event () {
  if (_time != TIMESTEP_NEVER)
    Simulation::gMainSummary->report(*this);
}


ostream& operator<<(ostream& out, const Event& event){
  out << event._time << endl;
  if (event._time == TIMESTEP_NEVER) return out;
  out << event._surveyPeriod << endl;
  out << event._ageGroup << endl;
  out << event._state << endl;
  out << event._outcome << endl;
  return out;
}

istream& operator>>(istream& in, Event& event){
  in >> event._time;
  if (event._time == TIMESTEP_NEVER) return in;
  in >> event._surveyPeriod;
  in >> event._ageGroup;
  int temp;
  in >> temp;
  event._state = Pathogenesis::State(temp);
  in >> event._outcome;
  return in;
}
