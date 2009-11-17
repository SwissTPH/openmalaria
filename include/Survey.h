/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef Hmod_Survey
#define Hmod_Survey

#include "Global.h"

/** Enumeration of reporting options
 *
 * Most are reported per age-group. */
enum Measure {
  /// Total number of humans
  nHost = 0,
  /// number of infected hosts 
  nInfect = 1,
  /// expected number of infected hosts
  nExpectd= 2,
  /// number of patent hosts
  nPatent= 3,
  /// Sum of the log of the pyrogen threshold
  sumLogPyrogenThres = 4,
  /// Sum of the logarithm of the parasite density
  sumlogDens= 5,
  /// Total infections
  totalInfs= 6,
  /** Infectiousness of human population to mosquitoes
   *
   * Number of hosts transmitting to mosquitoes (i.e. sum of proportion of
   * mosquitoes that get infected). We don't want this by age. */
  nTransmit= 7,
  /// Total patent infections
  totalPatentInf= 8,
  /// Contribution to immunity functions
  ///NOTE: not used
  contrib= 9,
  /// Sum of the pyrogenic threshold
  sumPyrogenThresh = 10,
  /// number of treatments (1st line)
  nTreatments1= 11,
  /// number of treatments (2nd line)
  nTreatments2= 12,
  /// number of treatments (inpatient)
  nTreatments3= 13,
  /// number of episodes (uncomplicated)
  nUncomp= 14,
  /// number of episodes (severe)
  nSevere= 15,
  /// cases with sequelae
  nSeq= 16,
  /// deaths in hospital
  nHospitalDeaths= 17,
  /// number of deaths (indirect)
  nIndDeaths= 18,
  /// number of deaths (direct)
  nDirDeaths= 19,
  /// number of EPI vaccine doses given
  nEPIVaccinations= 20,
  //all cause infant mortality rate
  imr_summary= 21,
  /// number of Mass / Campaign vaccine doses given
  nMassVaccinations= 22,
  /// recoveries in hospital
  nHospitalRecovs= 23,
  /// sequelae in hospital
  nHospitalSeqs= 24,
  /// number of IPT Doses
  nIPTDoses= 25,
  /** Annual Average Kappa
   *
   * Calculated once a year as sum of human infectiousness divided by initial
   * EIR summed over a year. */
  annAvgK= 26,
  /// Number of episodes (non-malaria fever)
  nNMFever= 27,
  /// Innoculations per human (all ages) per day of year, over the last year.
  innoculationsPerDayOfYear = 28,
  /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
  kappaPerDayOfYear = 29,
  /** The total number of innoculations per age group, summed over the
   * reporting period. */
  innoculationsPerAgeGroup = 30,
  // TODO: Replace encoding within XML with something more extensible.
  // Can't use 31 without being very careful, because we're using a *signed* 32-bit int.
  // Note: can't use values greater than 31 without forcing a 64-bit type
  NUM_SUMMARY_OPTIONS	// must be hightest value above plus one
};

/** Included for type-saftey: don't allow implicit double->int conversions.
 *
 * Incindentally, the constructor can be used implicitly for implicit
 * conversion doing the right thing.
 * 
 * Don't use _this_ class for other index/age-group types. */
class SurveyAgeGroup {
  public:
    /** Find the age group for the given age ageYears. */
    SurveyAgeGroup (double ageYears);
    /** Used for checkpointing. */
    void read (istream& in) {
      in >> _i;
    }
    
    /** Get the represented index. */
    inline size_t i () {
      return _i;
    }
    
    /// Get the total number of age categories (inc. one for indivs. not in any
    /// category given in XML).
    static inline int getNumGroups () {
      return _upperbound.size();
    }
    
  private:
    size_t _i;
    
    /// Initialize _lowerbound and _upperbound
    static void init ();
    
    static double _lowerbound; //!< Lower boundary of the youngest agegroup
    /** Upper boundary of agegroups, in years.
     *
     * These are age-groups given in XML plus one with no upper limit for
     * individuals outside other bounds. */
    static vector<double> _upperbound; 
    
    friend class Survey;
};

/// Data struct for a single survey.
class Survey {
  ///@brief Static members (options from XML)
  //@{
  private:
    /// Initialize static parameters.
    static void init();
    
    /// Encoding of which summary options are active in XML is converted into
    /// this array for easier reading (and to make changing encoding within XML easier).
    static bool active[NUM_SUMMARY_OPTIONS];
    
    /** Assimilator mode
     *
     * If true, skip the first 3 columns of output to reduce file size. */
    static bool _assimilatorMode; 
  //@}
  
public:
  /// @brief reportXXX functions to report val more of measure XXX within age-group ageGroup. Returns this allowing chain calling.
  //Note: generate this list from variable definitions by regexp search-replacing using the following:
  //Search: vector<(\w+)> _num(\w+)\;
  //Replace: Survey& report\2 (SurveyAgeGroup ageGroup, \1 val) {\n    _num\2[ageGroup.i()] += val;\n    return *this;\n  }
  //Search: vector<(\w+)> _sum(\w+)\;
  //Replace: Survey& addTo\2 (SurveyAgeGroup ageGroup, \1 val) {\n    _sum\2[ageGroup.i()] += val;\n    return *this;\n  }
  //@{
    Survey& reportHosts (SurveyAgeGroup ageGroup, int val) {
      _numHosts[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportInfectedHosts (SurveyAgeGroup ageGroup, int val) {
      _numInfectedHosts[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportExpectedInfected (SurveyAgeGroup ageGroup, double val) {
      _numExpectedInfected[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportPatentHosts (SurveyAgeGroup ageGroup, int val) {
      _numPatentHosts[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToLogPyrogenicThreshold (SurveyAgeGroup ageGroup, double val) {
      _sumLogPyrogenicThreshold[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToLogDensity (SurveyAgeGroup ageGroup, double val) {
      _sumLogDensity[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToInfections (SurveyAgeGroup ageGroup, int val) {
      _sumInfections[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToPatentInfections (SurveyAgeGroup ageGroup, int val) {
      _sumPatentInfections[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToPyrogenicThreshold (SurveyAgeGroup ageGroup, double val) {
      _sumPyrogenicThreshold[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportTreatments1 (SurveyAgeGroup ageGroup, int val) {
      _numTreatments1[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportTreatments2 (SurveyAgeGroup ageGroup, int val) {
      _numTreatments2[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportTreatments3 (SurveyAgeGroup ageGroup, int val) {
      _numTreatments3[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportUncomplicatedEpisodes (SurveyAgeGroup ageGroup, int val) {
      _numUncomplicatedEpisodes[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportSevereEpisodes (SurveyAgeGroup ageGroup, int val) {
      _numSevereEpisodes[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportSequelae (SurveyAgeGroup ageGroup, int val) {
      _numSequelae[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportHospitalDeaths (SurveyAgeGroup ageGroup, int val) {
      _numHospitalDeaths[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportIndirectDeaths (SurveyAgeGroup ageGroup, int val) {
      _numIndirectDeaths[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportDirectDeaths (SurveyAgeGroup ageGroup, int val) {
      _numDirectDeaths[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportEPIVaccinations (SurveyAgeGroup ageGroup, int val) {
      _numEPIVaccinations[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportMassVaccinations (SurveyAgeGroup ageGroup, int val) {
      _numMassVaccinations[ageGroup.i()] += val;
      return *this;
    } 
    Survey& reportHospitalRecoveries (SurveyAgeGroup ageGroup, int val) {
      _numHospitalRecoveries[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportHospitalSequelae (SurveyAgeGroup ageGroup, int val) {
      _numHospitalSequelae[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportIPTDoses (SurveyAgeGroup ageGroup, int val) {
      _numIPTDoses[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportNonMalariaFevers (SurveyAgeGroup ageGroup, int val) {
      _numNonMalariaFevers[ageGroup.i()] += val;
      return *this;
    } 
  //@}
  /// Report a first or second line, or inpatient treatment
  void reportTreatment(SurveyAgeGroup ageGroup, int regimen);
  
  void setAnnualAverageKappa(double kappa) {
    _annualAverageKappa = kappa;
  }
  void setNumTransmittingHosts(double value) {
    _numTransmittingHosts = value;
  }
  
  void setInnoculationsPerDayOfYear (vector<double>& v) {
    _innoculationsPerDayOfYear = v;
  }
  void setKappaPerDayOfYear (vector<double>& v) {
    _kappaPerDayOfYear = v;
  }
  void setInnoculationsPerAgeGroup (vector<double>& v) {
    _innoculationsPerAgeGroup = v;	// copies v, not just its reference
  }
  
private:
  /// Resize all vectors
  void allocate ();
  
  /** Write out arrays
   * @param outputFile Stream to write to
   * @param survey Survey number (starting from 1) */
  void writeSummaryArrays (ostream& outputFile, int survey);
  
  vector<int> _numHosts;
  vector<int> _numInfectedHosts;
  vector<double> _numExpectedInfected;
  vector<int> _numPatentHosts;
  vector<double> _sumLogPyrogenicThreshold;
  vector<double> _sumLogDensity;
  vector<int> _sumInfections;
  double _numTransmittingHosts;
  vector<int> _sumPatentInfections;
  vector<double> _sumPyrogenicThreshold;
  vector<int> _numTreatments1;
  vector<int> _numTreatments2;
  vector<int> _numTreatments3;
  vector<int> _numUncomplicatedEpisodes;
  vector<int> _numSevereEpisodes;
  vector<int> _numSequelae;
  vector<int> _numHospitalDeaths;
  vector<int> _numIndirectDeaths;
  vector<int> _numDirectDeaths;
  vector<int> _numEPIVaccinations;
  vector<int> _numMassVaccinations; 
  vector<int> _numHospitalRecoveries;
  vector<int> _numHospitalSequelae;
  vector<int> _numIPTDoses;
  double _annualAverageKappa;
  vector<int> _numNonMalariaFevers; 
  vector<double> _innoculationsPerDayOfYear;
  vector<double> _kappaPerDayOfYear;
  vector<double> _innoculationsPerAgeGroup;
  
  friend class SurveysType;
};

/** Line end character. Use Unix line endings to save a little size. */
const char lineEnd = '\n';

template <class T>
void writeArray(ostream& file, int measure, bool assimilationMode, int survey, vector<T>& array);

template <class T>
void writeArray(ostream& file, int measure, bool assimilationMode, int survey, T& value);

#endif
