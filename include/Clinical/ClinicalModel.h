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
#include "Episode.h"

namespace OM { namespace Clinical {

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
  
  /// Static checkpointing
  static void staticCheckpoint (istream& stream);
  static void staticCheckpoint (ostream& stream);
  
  /** Return a new ClinicalModel.
   *
   * @param cF 	comorbidity factor, passed to PathogenesisModel
   * @param tSF	treatment seeking factor, passed to CaseManagementModel */
  static ClinicalModel* createClinicalModel (double cF, double tSF);
  
  static void initMainSimulation ();
  
  /** Calculate infant mortality as deaths/1000 livebirths for the whole main-
   * simulation period (not as deaths/1000 years-at-risk per survey).
   * 
   * This mimicks field data on all-cause mortality in infants.
   * Uses the kaplan-meier method because the demography was set up to provide
   * a stable age-distribution but unfortunately does not accurately describe
   * death rates. The kaplan-meier estimate is the product of the proportion of
   * infants survivng at each interval. */
  static double infantAllCauseMort();
  //@}
  
  /// Destructor
  virtual ~ClinicalModel ();
  
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
  void update (WithinHost::WithinHostModel& withinHostModel, double ageYears, int ageTimeSteps);
  
  /** For infants, updates the infantIntervalsAtRisk and potentially
   * infantDeaths arrays. */
  void updateInfantDeaths (int ageTimeSteps);
  
  /** Was the last diagnosis severe malaria?
   * Used when mass treament is used with the DescriptiveIPT code. */
  inline bool latestDiagnosisIsSevereMalaria () {
    return latestReport.getState() == Pathogenesis::STATE_SEVERE;
  }
  
  /** Used with IPT within host model to potentially skip summarizing.
   *
   * NOTE: Only used for IPT, which can only be used with OldCaseManagement.
   * In other cases, this just returns false. */
  virtual bool recentTreatment() {
    return false;
  }
  
  /// Summarize PathogenesisModel details
  void summarize (Survey& survey, SurveyAgeGroup ageGroup);
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      checkpoint (stream);
  }
  
private:
    ///@brief Infant death summaries (checkpointed).
    //@{
    static vector<int> infantDeaths;
    static vector<int> infantIntervalsAtRisk;
    //@}
    
    /// Non-malaria mortality in under 1year olds.
    /// Set by init ()
    static double _nonMalariaMortality;
  
protected:
  /// Constructor.
  ClinicalModel (double cF);
  /// Constructor, loading from a checkpoint.
  ClinicalModel (istream& in);
  
  /** Update for clinical model - new pathogenesis status, treatment, etc.
   *
   * @param withinHostModel = WithinHostModel of human.
   * @param ageYears = Age of human. */
  virtual void doClinicalUpdate (WithinHost::WithinHostModel& withinHostModel, double ageYears) =0;
  
  virtual void checkpoint (istream& stream);
  virtual void checkpoint (ostream& stream);
  
  
  /// The PathogenesisModel introduces illness dependant on parasite density
  Pathogenesis::PathogenesisModel *pathogenesisModel;
  
  /** Next event to report.
   * Only reported when the Human dies or a separate episode occurs. */
  Episode latestReport;
  
  /** @brief Positive values of _doomed variable (codes). */
  enum {
    DOOMED_TOO_OLD = 1,		///< died because reached age limit
    DOOMED_COMPLICATED = 4,	///< died from severe malaria or malaria with a coinfection
    DOOMED_NEONATAL = 6,	///< died due to mother's malaria infection
    DOOMED_INDIRECT = 7,	///< died indirectly from malaria (after a delay)
  };
  
  /** Can indicate that the individual is dead or about to die.
   *
   * If _doomed < 0, the individual is doomed to die.
   * 
   * If _doomed > 0, the individual is dead, and will be removed from the
   * population at the beginning of the next timestep. NOTE: why not
   * immediately? See above enum for positive values used. */
  int _doomed;
};

} }
#endif
