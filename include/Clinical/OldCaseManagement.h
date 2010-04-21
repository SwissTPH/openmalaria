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
#ifndef Hcase_management
#define Hcase_management

#include "Clinical/CaseManagementCommon.h"
#include "WithinHost/WithinHostModel.h"
#include "Clinical/Episode.h"

namespace scnXml {
    class HSImmediateOutcomes;
}

namespace OM { namespace Clinical {

namespace Regimen {
    /** Regimen: UC / UC2 / SEVERE.
     *
     * Note: values used in array lookups, so are important. */
    enum Type {
	UC = 0,		// first line
	UC2 = 1,		// second line
	SEVERE = 2,	// third line
	NUM = 3,
    };
}

//! Models of treatment seeking and referral
class OldCaseManagement : public CaseManagementCommon {
public:
  /// Initialize static parameters
  static void init();
  
  /** Load health system data from initial data or an intervention's data (both from XML).
   * (Re)loads all data affected by this healthSystem element.
   *
   * @param source Values have same meaning as healthSystemSource variable. */
  static void setHealthSystem (int source);
  
  
  //!Read caseManagement parameters from input file and allocate data structures.
  OldCaseManagement(double tSF);
  /// Load from checkpoint
  OldCaseManagement(istream& in);
  ~OldCaseManagement();
  
  /** Determine treatment for a human.
   * @param pgState Wellbeing of subject (well, severe malaria sickness, etc.)
   * @param withinHostModel WithinHostModel of human.
   * @param latestReport Reporting memory
   * @param ageYears Age of human.
   * @param ageGroup Survey age group of human.
   * @param doomed _doomed variable of Human; used to kill the human.
   *	Passing like this isn't ideal. */
  void doCaseManagement (Pathogenesis::State pgState, WithinHost::WithinHostModel& withinHostModel, Episode& latestReport, double ageYears, SurveyAgeGroup ageGroup, int& doomed);
  
  inline bool recentTreatment() {
    return (Global::simulationTime-_tLastTreatment >= 1 &&
	    Global::simulationTime-_tLastTreatment <= 4);
  }
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      _tLastTreatment & stream;
      _treatmentSeekingFactor & stream;
  }
  
private:
   /** Called when a non-severe/complicated malaria sickness occurs.
    *
    * @returns True in case of effective or partially effective treatment, false otherwise. */
   bool uncomplicatedEvent(Episode& latestReport, bool isMalaria, double ageYears, SurveyAgeGroup ageGroup);

   /** Called when a severe/complicated (with co-infection) malaria sickness occurs.
    *
    * @returns True in case of effective or partially effective treatment, false otherwise.
    * 
    * Note: sets doomed = 4 if patient dies. */
  bool severeMalaria(Episode& latestReport, double ageYears, SurveyAgeGroup ageGroup, int& doomed);
  
  /** Timestep of the last treatment (TIMESTEP_NEVER if never treated). */
  int _tLastTreatment;
  
  //! treatment seeking for heterogeneity
  double _treatmentSeekingFactor;
  
  /*! Calculate the case fatality rate in the community as a function of the
    hospital case fatality rate.*/
  static double getCommunityCaseFatalityRate(double caseFatalityRatio);
  
  /// Calculate _probGetsTreatment, _probParasitesCleared and _cureRate.
  static void setParasiteCaseParameters (const scnXml::HSImmediateOutcomes& healthSystem);
  
  //log odds ratio of case-fatality in community compared to hospital
  //set only by init()
  static double _oddsRatioThreshold;
  
  //BEGIN Static parameters, set by setHealthSystem()
  // These parameters are reset via a setHealthSystem call on checkpoint
  // load rather than checkpointed.
  
  /// Age bounds of probSequelae* parameters
  //@{
  static const size_t NUM_SEQUELAE_AGE_GROUPS = 2;
  static const int SEQUELAE_AGE_BOUND[NUM_SEQUELAE_AGE_GROUPS];
  //@}
  
  /** pSequelaeTreated is the probability that the patient has sequelae
   * conditional on hospital treatment for severe disease. */
  static double probSequelaeTreated[NUM_SEQUELAE_AGE_GROUPS];
  /** pSequelaeUntreated is the probability that the patient has sequelae
   * conditional if they don't receive hospital treatment for severe disease.
   */
  static double probSequelaeUntreated[NUM_SEQUELAE_AGE_GROUPS];

  static double probGetsTreatment[Regimen::NUM];
  static double probParasitesCleared[Regimen::NUM];
  static double cureRate[Regimen::NUM];
  //END
};

} }
#endif
