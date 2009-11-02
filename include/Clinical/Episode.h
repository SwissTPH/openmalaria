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

#ifndef Hmod_Episode
#define Hmod_Episode

#include "Global.h"
#include <ostream>

/** Summary of clinical events during a caseManagementMemory period, in one individual.
 *
 * Terminology:
 * An "event" is an instantaeneous alteration of state.
 * An "episode" is a clinical-view of sickness caused by a malaria infection.
 * There's no reason an "episode" can't span multiple infections and multiple
 * bouts of sickness and recovery (the most severe is reported). */
class Episode{
public:
  Episode() : _time(TIMESTEP_NEVER) {};
  ~Episode();

  friend ostream& operator<<(ostream& out, const Episode& event);
  friend istream& operator>>(istream& in, Episode& event);

  /** Report an episode, its severity, and any outcomes it entails.
   *
   * @param simulationTime Time of report (i.e. now)
   * @param ageGroup Monitoring agegroup
   * @param newState The severity (diagnosis) and outcome.
   */
  void update(int simulationTime, int ageGroup, Pathogenesis::State newState);

  Pathogenesis::State getState() const {return _state;};
  int getAgeGroup() const {return _ageGroup;};
  int getSurveyPeriod() const {return _surveyPeriod;};
  
private:
  /// Timestep of event (TIMESTEP_NEVER if no event).
  int _time;
  //! survey period during which the event occured
  //! TODO: we could use the survey array to map time to survey period. slower, but less memory.
  int _surveyPeriod;
  //! agegroup of the individual which experienced the episode
  int _ageGroup;
  /// Descriptor of state, containing reporting info. Not all information will
  /// be reported (e.g. indirect deaths are reported independantly).
  Pathogenesis::State _state;
};

#endif
