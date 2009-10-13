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
#ifndef Hmod_event
#define Hmod_event
#include <fcntl.h>
#include <sys/types.h>
#include <math.h>
#include "Global.h"
#include <ostream>

/*
  Summary for clinical events during a caseManagementMemory period, in one individual.
  time in 5 day intervals
*/

class Event{
public:
  Event() : _time(TIMESTEP_NEVER) {};
  ~Event();

  friend ostream& operator<<(ostream& out, const Event& event);
  friend istream& operator>>(istream& in, Event& event);

  //! Report and replace a previous episode, or update the diagnosis/outcome.
  /*!
    \param agegroup Monitoring agegroup
    \param diagnosis
    \outcome 
  */
  void update(int simulationTime, int ageGroup, int diagnosis, int outcome);

  int getDiagnosis() const {return _diagnosis;};
  int getAgeGroup() const {return _ageGroup;};
  int getSurveyPeriod() const {return _surveyPeriod;};
  int getOutcome() const {return _outcome;};
  
private:
  /// Timestep of event (TIMESTEP_NEVER if no event).
  int _time;
  //! survey period during which the event occured
  //! TODO: we could use the survey array to map time to survey period. slower, but less memory.
  int _surveyPeriod;
  //! agegroup of the individual which experienced the episode
  int _ageGroup;
  //! final diagnosis, severe if one of the clinical events was severe, else uncomplicated
  int _diagnosis;
  //! maximum of recovered, sequelae, death
  int _outcome;
  //! The total(!) number of clinical events that occured during this caseManagementMemory period
  int _recurrence;
};

#endif
