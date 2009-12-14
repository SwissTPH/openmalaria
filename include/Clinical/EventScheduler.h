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
#include <map>
#include <list>

namespace Clinical {

/** Tracks clinical status (sickness), triggers case management for new events,
 * medicates treatment, determines patient recovery, death and sequelae.
 * 
 * TODO: Reporting of parasitological status (not model specific).
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
  /// Current state of sickness
  Pathogenesis::State pgState;
  /// Time of last state-change; only meaningful if pgState & Pathogenesis::SICK.
  int pgChangeTimestep;
  
  /// All pending medications
  list<MedicateData> medicateQueue;
  /// Decision ID of last case management run
  cmid lastCmDecision;
  
  ///@brief Static data, set up by init
  //@{
  //FIXME: use hash-map with ID from DecisionEnums
  /*const static size_t PTABLE_SIZE = TREATMENT_NUM_TYPES * DRUG_NUM_SEV;
  static double pDeathTable[PTABLE_SIZE];
  static double pRecoverTable[PTABLE_SIZE];*/
  
  //@}
};

}
#endif
