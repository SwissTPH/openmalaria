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
#include "summary.h"

int Episode::reportingPeriodMemory;


void Episode::update (int simulationTime, int ageGroup, Pathogenesis::State newState)
{
  if (simulationTime > (_time + reportingPeriodMemory)) {
    if (_time != TIMESTEP_NEVER) {
      Simulation::gMainSummary->report (*this);
    }
    _time = simulationTime;
    _surveyPeriod = Simulation::gMainSummary->getSurveyPeriod();
    _ageGroup = ageGroup;
    _state = newState;
  } else {
    _state = Pathogenesis::State (_state | newState);
  }
}

Episode::~Episode ()
{
  if (_time != TIMESTEP_NEVER)
    Simulation::gMainSummary->report (*this);
}


ostream& operator<< (ostream& out, const Episode& event)
{
  out << event._time << endl;
  if (event._time == TIMESTEP_NEVER) return out;
  out << event._surveyPeriod << endl;
  out << event._ageGroup << endl;
  out << event._state << endl;
  return out;
}

istream& operator>> (istream& in, Episode& event)
{
  in >> event._time;
  if (event._time == TIMESTEP_NEVER) return in;
  in >> event._surveyPeriod;
  in >> event._ageGroup;
  int temp;
  in >> temp;
  event._state = Pathogenesis::State (temp);
  return in;
}
