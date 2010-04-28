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

#ifndef Hmod_ClinicalEventSchduler
#define Hmod_ClinicalEventSchduler

#include "Global.h"
#include "Clinical/ClinicalModel.h"
#include "Clinical/ESCaseManagement.h"
#include <boost/unordered_map.hpp>
#include <list>

namespace scnXml {
    class HSEventScheduler;
}

namespace OM { namespace Clinical {

/** Tracks clinical status (sickness), triggers case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 * 
 * TODO: Reporting of parasitological status (not model specific).
 * 
 * Note: there are several variables that only need to be used during an
 * episode. It's possible that memory usage could be reduced by storing them
 * externally in a temporary object during episodes (but unlikely worth doing).
 */
class ClinicalEventScheduler : public ClinicalModel
{
public:
  static void init ();
  static void setParameters (const scnXml::HSEventScheduler& esData);
  static void cleanup ();
  
  ClinicalEventScheduler (double cF, double tSF);
  ~ClinicalEventScheduler ();
  
  virtual void massDrugAdministration(WithinHost::WithinHostModel& withinHostModel, double ageYears);
  
protected:
  virtual void doClinicalUpdate (WithinHost::WithinHostModel& withinHostModel, PerHostTransmission& hostTransmission, double ageYears, SurveyAgeGroup ageGroup);
  
  virtual void checkpoint (istream& stream);
  virtual void checkpoint (ostream& stream);
 
private:
    /// Utility find to take a value from pDeathInitial.
    double getPDeathInitial (double ageYears);
    
    /// Maximum number of timesteps (including first of case) an individual will
    /// remember they are sick before resetting.
    static int maxUCSeekingMemory;
    /// Length of an uncomplicated case
    static int uncomplicatedCaseDuration;
    /// Length of a complicated case
    static int complicatedCaseDuration;
    /// Time-span for which individual is at risk of death in complicated case
    /// minus length of complicated case (must be <= 0)
    static int extraDaysAtRisk;
    /// Probability that UC treatment seeking will be done immediately when
    /// sick, on second day given that it wasn't done on first, etc.
    static double pImmediateUC;
    
    /// Probability of death on first day of a complicated case (1 - S(0)).
    /// Depends on age-group.
    static map<double,double> pDeathInitial;
    
    /// Parameter of S(t) for t > 0
    static double neg_v;
    
    /// Multiplies the odds of dying in a community case, where the base
    /// (hospital) probability of death comes from other inputs.
    static double communityOddsMultiplier;
    
    // Note on memory usage: Pathogenesis::State is and enum (an int), so we
    // have a vtable followed by 3 ints, a double and a list. Alignment probably
    // wastes some space.
  /// Current state of sickness
  Pathogenesis::State pgState;
  
  //NOTE: not all these variables should be needed eventually
  /** Set to when an event should start. If simulationTime equals this, an event
   * is started (UC & C behaviour different).
   * 
   * Note: medications are not delayed by this. */
  int caseStartTime;
  
  /** The individual recovers when Global::simulationTime >= timeOfRecovery,
   * assuming they didn't die. */
  int timeOfRecovery;
  
  /// Total parasite density at previous timestep (used during an event).
  double previousDensity;
  
  /// All pending medications
  list<MedicateData> medicateQueue;
};

} }
#endif