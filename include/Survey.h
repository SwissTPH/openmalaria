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
#include <bitset>
#include <map>

namespace OM {
    
/** Enumeration of reporting options
 *
 * Many are reported per age-group, but to check which actually are you'll have
 * to look through the code.
 * 
 * Don't ever change these names or numbers. The names are used in scenario
 * files, and the numbers in results output/databases. */
enum SurveyMeasure {
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
     * mosquitoes that get infected). Single value, not per age-group. */
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
     * EIR summed over a year. Single value, not per age-group. */
    annAvgK= 26,
    /// Number of episodes (non-malaria fever)
    nNMFever= 27,
    /** The total number of innoculations per age group, summed over the
     * reporting period. */
    innoculationsPerAgeGroup = 30,
    
    //BEGIN Per day-of-year data (removed)
    /// Innoculations per human (all ages) per day of year, over the last year.
    innoculationsPerDayOfYear = 28,
    /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
    kappaPerDayOfYear = 29,
    //END
    
    /** @brief Vector model parameters.
     *
     * All are point-time outputs, not averages. The Nv0, Nv, Ov and Sv outputs
     * are per-species; the EIR outputs are single values. */
    //@{
    /** Mosquito emergence rate. */
    Vector_Nv0 = 31,
    /// Mosquito population size
    Vector_Nv = 32,
    /// Number of infected mosquitoes
    Vector_Ov = 33,
    /// Number of infectious mosquitoes
    Vector_Sv = 34,
    /// Input EIR (Expected EIR entered into scenario file)
    Vector_EIR_Input = 35,
    /// Simulated EIR (EIR output by the vector model)
    Vector_EIR_Simulated = 36,
    //@}
    
    // must be hightest value above plus one
    NUM_SURVEY_OPTIONS	
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
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
	_i & stream;
    }
    
    /** Get the represented index. */
    inline size_t i () {
      return _i;
    }
    
    /// Get the total number of age categories (inc. one for indivs. not in any
    /// category given in XML).
    static inline size_t getNumGroups () {
      return _upperbound.size();
    }
    
  private:
    size_t _i;
    
    /// Initialize _lowerbound and _upperbound
    static void init ();
    
    //BEGIN Static parameters only set by init()
    /// Lower boundary of the youngest agegroup
    static double _lowerbound;
    /** Upper boundary of agegroups, in years.
     *
     * These are age-groups given in XML plus one with no upper limit for
     * individuals outside other bounds. */
    static vector<double> _upperbound;
    //END
    
    friend class Survey;
};

/// Data struct for a single survey.
class Survey {
  ///@brief Static members (options from XML). Parameters only set by init().
  //@{
  private:
    /// Initialize static parameters.
    static void init();
    
    /// Encoding of which summary options are active in XML is converted into
    /// this array for easier reading (and to make changing encoding within XML easier).
    static bitset<NUM_SURVEY_OPTIONS> active;
    
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
  
  void setInnoculationsPerAgeGroup (vector<double>& v) {
    _innoculationsPerAgeGroup = v;	// copies v, not just its reference
  }
  void reportHospitalizationDays (int days) {
      //FIXME: report
  }
  void reportHospitalEntries (int num) {
      //FIXME: report
  }
  void reportRDT (int num) {
      //FIXME: report
  }
  void set_Vector_Nv0 (string key, double v) {
    data_Vector_Nv0[key] = v;
  }
  void set_Vector_Nv (string key, double v) {
    data_Vector_Nv[key] = v;
  }
  void set_Vector_Ov (string key, double v) {
    data_Vector_Ov[key] = v;
  }
  void set_Vector_Sv (string key, double v) {
    data_Vector_Sv[key] = v;
  }
  void set_Vector_EIR_Input (double v) {
    data_Vector_EIR_Input = v;
  }
  void set_Vector_EIR_Simulated (double v) {
    data_Vector_EIR_Simulated = v;
  }
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
    _numHosts & stream;
    _numInfectedHosts & stream;
    _numExpectedInfected & stream;
    _numPatentHosts & stream;
    _sumLogPyrogenicThreshold & stream;
    _sumLogDensity & stream;
    _sumInfections & stream;
    _numTransmittingHosts & stream;
    _sumPatentInfections & stream;
    _sumPyrogenicThreshold & stream;
    _numTreatments1 & stream;
    _numTreatments2 & stream;
    _numTreatments3 & stream;
    _numUncomplicatedEpisodes & stream;
    _numSevereEpisodes & stream;
    _numSequelae & stream;
    _numHospitalDeaths & stream;
    _numIndirectDeaths & stream;
    _numDirectDeaths & stream;
    _numEPIVaccinations & stream;
    _numMassVaccinations & stream; 
    _numHospitalRecoveries & stream;
    _numHospitalSequelae & stream;
    _numIPTDoses & stream;
    _annualAverageKappa & stream;
    _numNonMalariaFevers & stream; 
    _innoculationsPerAgeGroup & stream;
    data_Vector_Nv0 & stream;
    data_Vector_Nv & stream;
    data_Vector_Ov & stream;
    data_Vector_Sv & stream;
    data_Vector_EIR_Input & stream;
    data_Vector_EIR_Simulated & stream;
  }
  
private:
  /// Resize all vectors
  void allocate ();
  
  /** Write out arrays
   * @param outputFile Stream to write to
   * @param survey Survey number (starting from 1) */
  void writeSummaryArrays (ostream& outputFile, int survey);
  
  // atomic data:
  double _numTransmittingHosts;
  double _annualAverageKappa;
  
  // data, per SurveyAgeGroup:
  vector<int> _numHosts;
  vector<int> _numInfectedHosts;
  vector<double> _numExpectedInfected;
  vector<int> _numPatentHosts;
  vector<double> _sumLogPyrogenicThreshold;
  vector<double> _sumLogDensity;
  vector<int> _sumInfections;
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
  vector<int> _numNonMalariaFevers; 
  vector<double> _innoculationsPerAgeGroup;
  
    // data, per vector species:
    map<string,double> data_Vector_Nv0;
    map<string,double> data_Vector_Nv;
    map<string,double> data_Vector_Ov;
    map<string,double> data_Vector_Sv;
    double data_Vector_EIR_Input;
    double data_Vector_EIR_Simulated;
  
  friend class SurveysType;
};

/** Line end character. Use Unix line endings to save a little size. */
const char lineEnd = '\n';

}
#endif
