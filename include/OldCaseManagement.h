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

#include "CaseManagementModel.h"

//! Models of treatment seeking and referral
class OldCaseManagement : public CaseManagementModel {
public:
  /// Initialize static parameters
  static void init();
  
  //!Read caseManagement parameters from input file and allocate data structures.
  OldCaseManagement(double tSF);
  /// Load from checkpoint
  OldCaseManagement(istream& in);
  virtual ~OldCaseManagement();
  
  virtual void doCaseManagement (Pathogenesis::Infection, WithinHostModel&, Event& latestReport, double, int& doomed);
  
  virtual  void write(ostream& out) const;
  
private:
   /*! should return true in case of effective or partially effective
  treatment, false otherwise */
   bool uncomplicatedEvent(Event& latestReport, bool isMalaria, double ageYears);

  //! returns true in case of effective treatment, false otherwise
  bool severeMalaria(Event& latestReport, double ageYears, int& doomed);
  
  //!indicates the latest treatment regimen(1st, 2nd or 3rd line)
  int _latestRegimen;
  
  /*! Linear interpolation to get age-specific hospital case fatality rates
   * 
   * @param ageyears Age of person in years */
  static double caseFatality(double ageyears);

  /*! Calculate the case fatality rate in the community as a function of the
    hospital case fatality rate.*/
  static double getCommunityCaseFatalityRate(double caseFatalityRatio);

  /*!
    Look up any recent treatments and determine which drug regimen to use
    next.  latestTreatment: time of the most recent treatment for the
    individual regimen: drug to be used for the new treatment.
  */
  static int getNextRegimen(int simulationTime, int diagnosis, int tLastTreated, int regimen);

  /// Calculate _probGetsTreatment, _probParasitesCleared and _cureRate.
  //@{
  static void setParasiteCaseParameters ();
  //@}
  
  static double probGetsTreatment[3];
  static double probParasitesCleared[3];
  static double cureRate[3];
  
  //log odds ratio of case-fatality in community compared to hospital
  static double _oddsRatioThreshold;
  
  /// Age bounds of probSequelae* parameters
  //@{
  static const int NUM_SEQUELAE_AGE_GROUPS = 2;
  static const int SEQUELAE_AGE_BOUND[NUM_SEQUELAE_AGE_GROUPS];
  //@}
  
  /** pSequelaeTreated is the probability that the patient has sequelae
   * conditional on hospital treatment for severe disease. */
  static double probSequelaeTreated[2];
  /** pSequelaeUntreated is the probability that the patient has sequelae
   * conditional if they don't receive hospital treatment for severe disease.
   */
  static double probSequelaeUntreated[2];

  /// shortcut: if there is only one CFR group, and the CFR is 0, set this to true.
  static bool _noMortality;
  
  /// Age-specific bounds and case-fatality rates.
  //@{
  /// Age groups have the bounds [_inputAge[i], _inputAge[i+1])
  static std::vector<double> _inputAge;
  /** Case fatality rate for age groups; last entry is a copy of the previous
   * entry. */
  static std::vector<double> _caseFatalityRate;
  //@}

  //! Reads in the Case Fatality percentages from the XML.
  /*! This replaces the reading from CaseFatalityByAge.csv.Note that we could
    calculate and cache the CFR as a function of age in years for better
    performance. This would require a specification of the resolution.
  */
  static void readCaseFatalityRatio();
};

#endif
