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

#include <vector>

//! Models of treatment seeking and referral
class CaseManagement{

 public:

  //!Read caseManagement parameters from input file and allocate data structures.
  CaseManagement();
  ~CaseManagement();

  /*!
    Linear interpolation to get age-specific hospital case fatality rates
    ageyears: age of person in years
  */
  double caseFatality(double ageyears);

  /*! Calculate the case fatality rate in the community as a function of the
    hospital case fatality rate.*/
  double getCommunityCaseFatalityRate(double caseFatalityRatio) const;

  /*!
    Look up any recent treatments and determine which drug regimen to use
    next.  latestTreatment: time of the most recent treatment for the
    individual regimen: drug to be used for the new treatment.
  */
  int getNextRegimen(int simulationTime, int diagnosis, int tLastTreated, int regimen) const;
  bool indirectDeath(int dob, int ageGroup, int& doomed);

  int getCaseManagementMemory() const;
  double getRiskFromMaternalInfection() const;
  void setRiskFromMaternalInfection(double risk);

  double getProbabilityGetsTreatment(int regimen) const;
  double getProbabilityParasitesCleared(int regimen) const;
  double getCureRate(int regimen) const;
  double getProbabilitySequelaeTreated(int regimen) const;
  double getProbabilitySequelaeUntreated(int regimen) const;

 private:

  double _probGetsTreatment[3];
  double _probParasitesCleared[3];
  double _cureRate[3];
  int _caseManagementMemory;
  /*
    Probability for a newborn to die (indirect death) because the mother is infected.
    Depends on the prevalence of parasitaemia in mother at some previous t.
  */
  double _riskFromMaternalInfection;
  //log odds ratio of case-fatality in community compared to hospital
  double _oddsRatioThreshold;
  /*
    pSequelaeTreated is the probability that the patient has sequelae conditional on hospital
    treatment for severe disease. pSequelaeUntreated is the probability that the patient has 
    sequelae conditional if they don^t receive hospital treatment for severe disease.
  */
  double _probSequelaeTreated[2];
  double _probSequelaeUntreated[2];
  //shortcut logical: if there is only one CFR group, and the CFR is 0, set this to .true.
  bool _noMortality;
  //The next 2 arrays contain the age-specific case-fatality rates.
  std::vector<double> _inputAge;
  std::vector<double> _caseFatalityRate;

  //! Reads in the Case Fatality percentages from the XML.
  /*! This replaces the reading from CaseFatalityByAge.csv.Note that we could
    calculate and cache the CFR as a function of age in years for better
    performance. This would require a specification of the resolution.
  */
  void readCaseFatalityRatio();

};

#endif
