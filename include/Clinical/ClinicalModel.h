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

#ifndef Hmod_ClinicalModel
#define Hmod_ClinicalModel

#include "Pathogenesis/PathogenesisModel.h"
#include "event.h"

/** The clinical model models the effects of sickness dependant on malarial
 * parasite densities and administers anti-malaria treatments via the drug
 * model (or in a simpler case, directly clearing infections).
 * 
 * So far, sickness types include uncomplicated and severe malaria cases and
 * non-malaria sickness.
 * 
 * Patient outcomes include full recovery, recovery with sequelae and death.
 * 
 * Reporting includes patient outcome and potentially drug usage and use of
 * RDTs (Rapid Diagnostic Tests) for costing purposes. */
class ClinicalModel
{
public:
  /// @brief Static functions
  //@{
  /// Initialise whichever model is in use.
  static void init ();
  
  /** Return a new ClinicalModel.
   *
   * @param cF 	comorbidity factor, passed to PathogenesisModel
   * @param tSF	treatment seeking factor, passed to CaseManagementModel */
  static ClinicalModel* createClinicalModel (double cF, double tSF);
  /** Load a ClinicalModel from a checkpoint. */
  static ClinicalModel* createClinicalModel (istream& in);
  //@}
  
  /// Destructor
  virtual ~ClinicalModel ();
  /// Write a checkpoint
  virtual void write (ostream& out) =0;
  
  /** Kills the human if ageTimeSteps reaches the simulation age limit.
   *
   * @returns True if the human has been killed by some means. The clinical
   *	model now tracks this status. */
  bool isDead (int ageTimeSteps);
  
  /** Run main part of the model: determine the sickness status and any
   * treatment for the human.
   * 
   * @param withinHostModel = Used to get the parasite density and to medicate
   *	drugs/clear infections.
   * @param ageYears = Age of human. */
  void update (WithinHostModel& withinHostModel, double ageYears, int ageTimeSteps);
  
  /** For infants, updates the Global::infantIntervalsAtRisk and potentially
   * Global::infantDeaths arrays. */
  void updateInfantDeaths (int ageTimeSteps);
  
  /** Was the last diagnosis severe malaria?
   * FIXME: update or remove */
  inline bool latestDiagnosisIsSevereMalaria () {
    return latestReport.getDiagnosis() == Diagnosis::SEVERE_MALARIA;
  }
  
  /** Used with IPT within host model to potentially skip summarizing.
   *
   * NOTE: Only used for IPT, which can only be used with OldCaseManagement.
   * In other cases, this just returns false. */
  virtual bool recentTreatment() {
    return false;
  }
  
  /// Summarize PathogenesisModel details
  void summarize (Summary& summary, double age);
  
  /** The maximum age, in timesteps, of when a sickness event occurred, for
   * another event to be considered part of the same reported "event".
   * 
   * Used by both the clinical models in roughly the same way, but will have
   * different values in each to match Global::interval.
   * 
   * NOTE: notation: episode/event */
  static int reportingPeriodMemory;
  
protected:
  /// Constructor.
  ClinicalModel (double cF);
  /// Constructor, loading from a checkpoint.
  ClinicalModel (istream& in);
  
  /** Update for clinical model - new pathogenesis status, treatment, etc.
   *
   * @param withinHostModel = WithinHostModel of human.
   * @param ageYears = Age of human. */
  virtual void doClinicalUpdate (WithinHostModel& withinHostModel, double ageYears) =0;
  
  /// The PathogenesisModel introduces illness dependant on parasite density
  PathogenesisModel *pathogenesisModel;
  
  /** Next event to report.
   * Only reported when the Human dies or a separate episode occurs. */
  Event latestReport;
  
  /** Can indicate that the individual is dead or about to die.
   *
   * If _doomed < 0, the individual is doomed to die.
   * 
   * If _doomed > 0, the individual is dead, and will be removed from the
   * population at the beginning of the next timestep. NOTE: why not
   * immediately? Several positive values are used in the code. FIXME: why
   * is _doomed tested for values 4, -30 (5 days before indirect death) and
   * 6?
   * 4 = Direct death from severe malaria or a coinfection
   * 6 = Neonatal death
   * 7 = Indirect death */
  int _doomed;
};
#endif
