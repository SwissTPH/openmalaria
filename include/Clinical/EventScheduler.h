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
#include "Clinical/DecisionEnums.d"
#include <map>
#include <list>

namespace scnXml {
  class CaseType;
}

/** Tracks clinical status (sickness), does case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 * 
 * TODO: Case management cleanup
 * Case management delayed calls to medicate(), to allow cancelling when
 * patient becomes severe (from uncomplicated) or dies.
 * Reporting of parasitological status (not model specific).
 */
class ClinicalEventScheduler : public ClinicalModel
{
public:
  static void init ();
  
  ClinicalEventScheduler (double cF, double tSF);
  ClinicalEventScheduler (istream& in);
  ~ClinicalEventScheduler ();
  
  void write (ostream& out);
  
  void doClinicalUpdate (WithinHostModel& withinHostModel, double ageYears);
  
private:
  void doCaseManagement (WithinHostModel& withinHostModel, double ageYears);
  
  /// Current state of sickness
  Pathogenesis::State pgState;
  /// Time of last state-change; only meaningful if pgState & Pathogenesis::SICK.
  int pgChangeTimestep;
  
  /// Data used for a withinHostModel->medicate() call
  struct MedicateData {
    string abbrev;	/// Drug abbreviation
    double qty;		/// Quantity of drug prescribed
    int time;		/// Time of day to medicate at (minutes from start)
    int seekingDelay;	/// Delay before treatment seeking in days
  };
  
  /// All pending medications
  list<MedicateData> medicateQueue;
  /// Decision ID of last case management run
  int lastCmDecision;
  
  /** Cumulative probabilities and decisions for a case type */
  struct CaseTypeEndPoints {
    vector<double> cumProbs;
    vector<int> decisions;
  };
  
  /// Data type stored in decisions
  /// TODO: treatment seeking delay(?), hospital/community care, RDTs or not.
  struct CaseTreatment {
    /// Data for each medicate() call.
    vector<MedicateData> medications;
  };
  
  struct CaseManagementEndPoints {
    CaseTypeEndPoints caseUC1;
    CaseTypeEndPoints caseUC2;
    CaseTypeEndPoints caseSev;
    CaseTypeEndPoints caseNMFWithParasites;
    CaseTypeEndPoints caseNMFWithoutParasites;
    
    /** Decisions data (end points of decision tree plus part of path).
     *
     * NOTE: Using a red-black tree here; probably a hash-map would be faster.
     * NOTE: Decisions from all age groups are being combined. */
    map<size_t,CaseTreatment> decisions;
  };
  
  const static size_t PTABLE_NUM_DAYS = 3;
  static double pDeathTable[TREATMENT_NUM_TYPES * PTABLE_NUM_DAYS];
  static double pRecoverTable[TREATMENT_NUM_TYPES * PTABLE_NUM_DAYS];
  
  /** Age groups */
  static vector<double> caseManagementMaxAge;
  /** Case management data, per age group. */
  static vector<CaseManagementEndPoints> caseManagementEndPoints;
  
  static CaseTypeEndPoints readEndPoints (const scnXml::CaseType& caseType);
};
#endif
