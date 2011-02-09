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

#include "Monitoring/SurveyMeasure.h"
#include "Global.h"
#include <bitset>
#include <map>

namespace OM { namespace Monitoring {

/** Included for type-saftey: don't allow implicit double->int conversions.
 *
 * Incindentally, the constructor can be used implicitly for implicit
 * conversion doing the right thing.
 * 
 * Don't use _this_ class for other index/age-group types. */
class AgeGroup {
  public:
    AgeGroup () : _i(0) {}
    
    /** Update age-group. Assumes age only increases (per instance).
     *
     * If called regularly, should be O(1); worst case is O(_upperbound.size()). */
    void update (double ageYears);
    
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
    static bitset<SM::NUM_SURVEY_OPTIONS> active;
  //@}
  
public:
  /// @brief reportXXX functions to report val more of measure XXX within age-group ageGroup. Returns this allowing chain calling.
  //Note: generate this list from variable definitions by regexp search-replacing using the following:
  //Search: vector<(\w+)> _num(\w+)\;
  //Replace: Survey& report\2 (AgeGroup ageGroup, \1 val) {\n        _num\2[ageGroup.i()] += val;\n        return *this;\n    }
  //Search: vector<(\w+)> _sum(\w+)\;
  //Replace: Survey& addTo\2 (AgeGroup ageGroup, \1 val) {\n        _sum\2[ageGroup.i()] += val;\n        return *this;\n    }
  //@{
    Survey& reportHosts (AgeGroup ageGroup, int val) {
      _numHosts[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportInfectedHosts (AgeGroup ageGroup, int val) {
      _numInfectedHosts[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportExpectedInfected (AgeGroup ageGroup, double val) {
      _numExpectedInfected[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportPatentHosts (AgeGroup ageGroup, int val) {
      _numPatentHosts[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToLogPyrogenicThreshold (AgeGroup ageGroup, double val) {
      _sumLogPyrogenicThreshold[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToLogDensity (AgeGroup ageGroup, double val) {
      _sumLogDensity[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToInfections (AgeGroup ageGroup, int val) {
      _sumInfections[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToPatentInfections (AgeGroup ageGroup, int val) {
      _sumPatentInfections[ageGroup.i()] += val;
      return *this;
    }
    Survey& addToPyrogenicThreshold (AgeGroup ageGroup, double val) {
      _sumPyrogenicThreshold[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportTreatments1 (AgeGroup ageGroup, int val) {
      _numTreatments1[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportTreatments2 (AgeGroup ageGroup, int val) {
      _numTreatments2[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportTreatments3 (AgeGroup ageGroup, int val) {
      _numTreatments3[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportUncomplicatedEpisodes (AgeGroup ageGroup, int val) {
      _numUncomplicatedEpisodes[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportSevereEpisodes (AgeGroup ageGroup, int val) {
      _numSevereEpisodes[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportSequelae (AgeGroup ageGroup, int val) {
      _numSequelae[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportHospitalDeaths (AgeGroup ageGroup, int val) {
      _numHospitalDeaths[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportIndirectDeaths (AgeGroup ageGroup, int val) {
      _numIndirectDeaths[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportDirectDeaths (AgeGroup ageGroup, int val) {
      _numDirectDeaths[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportEPIVaccinations (AgeGroup ageGroup, int val) {
      _numEPIVaccinations[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportMassVaccinations (AgeGroup ageGroup, int val) {
      _numMassVaccinations[ageGroup.i()] += val;
      return *this;
    } 
    Survey& reportHospitalRecoveries (AgeGroup ageGroup, int val) {
      _numHospitalRecoveries[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportHospitalSequelae (AgeGroup ageGroup, int val) {
      _numHospitalSequelae[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportIPTDoses (AgeGroup ageGroup, int val) {
      _numIPTDoses[ageGroup.i()] += val;
      return *this;
    }
    Survey& reportNonMalariaFevers (AgeGroup ageGroup, int val) {
      _numNonMalariaFevers[ageGroup.i()] += val;
      return *this;
    } 
    Survey& reportNewInfections (AgeGroup ageGroup, int val) {
	_numNewInfections[ageGroup.i()] += val;
	return *this;
    }
    Survey& reportMassITNs (AgeGroup ageGroup, int val) {
	_numMassITNs[ageGroup.i()] += val;
	return *this;
    }
    Survey& reportEPI_ITNs (AgeGroup ageGroup, int val) {
	_numEPI_ITNs[ageGroup.i()] += val;
	return *this;
    }
    Survey& reportMassIRS (AgeGroup ageGroup, int val) {
	_numMassIRS[ageGroup.i()] += val;
	return *this;
    }
    Survey& reportMassVA (AgeGroup ageGroup, int val) {
	_numMassVA[ageGroup.i()] += val;
	return *this;
    }
    Survey& reportAddedToCohort (AgeGroup ageGroup, int val) {
        _numAddedToCohort[ageGroup.i()] += val;
        return *this;
    }
    Survey& reportRemovedFromCohort (AgeGroup ageGroup, int val) {
        _numRemovedFromCohort[ageGroup.i()] += val;
        return *this;
    }
    Survey& reportMDA (AgeGroup ageGroup, int val) {
        _numMDAs[ageGroup.i()] += val;
        return *this;
    }
    Survey& reportNmfDeaths (AgeGroup ageGroup, int val) {
        _numNmfDeaths[ageGroup.i()] += val;
        return *this;
    }
    Survey& reportAntibioticTreatments (AgeGroup ageGroup, int val) {
        _numAntibioticTreatments[ageGroup.i()] += val;
        return *this;
    }
  //@}
  
  void setAnnualAverageKappa(double kappa) {
    _annualAverageKappa = kappa;
  }
  void setNumTransmittingHosts(double value) {
    _numTransmittingHosts = value;
  }
  
  void setInoculationsPerAgeGroup (vector<double>& v) {
    _inoculationsPerAgeGroup = v;	// copies v, not just its reference
  }
  void report_Clinical_RDTs (int num) {
      _numClinical_RDTs += num;
  }
  void report_Clinical_DrugUsage (string abbrev, double qty) {
      // Insert the pair (abbrev, 0.0) if not there, get an iterator to it, and increment it's second param (quantity) by qty
      (*((_sumClinical_DrugUsage.insert(make_pair(abbrev, 0.0))).first)).second += qty;
  }
  void report_Clinical_DrugUsageIV (string abbrev, double qty) {
      // Insert the pair (abbrev, 0.0) if not there, get an iterator to it, and increment it's second param (quantity) by qty
      (*((_sumClinical_DrugUsageIV.insert(make_pair(abbrev, 0.0))).first)).second += qty;
  }
  Survey& report_Clinical_FirstDayDeaths (AgeGroup ageGroup, int val) {
      _numClinical_FirstDayDeaths[ageGroup.i()] += val;
      return *this;
  } 
  Survey& report_Clinical_HospitalFirstDayDeaths (AgeGroup ageGroup, int val) {
      _numClinical_HospitalFirstDayDeaths[ageGroup.i()] += val;
      return *this;
  } 
  void report_Clinical_Microscopy (int num) {
      _numClinical_Microscopy += num;
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
    _inoculationsPerAgeGroup & stream;
    data_Vector_Nv0 & stream;
    data_Vector_Nv & stream;
    data_Vector_Ov & stream;
    data_Vector_Sv & stream;
    data_Vector_EIR_Input & stream;
    data_Vector_EIR_Simulated & stream;
    _numClinical_RDTs & stream;
    _sumClinical_DrugUsage & stream;
    _sumClinical_DrugUsageIV & stream;
    _numClinical_FirstDayDeaths & stream;
    _numClinical_HospitalFirstDayDeaths & stream;
    _numNewInfections & stream;
    _numMassITNs & stream;
    _numEPI_ITNs & stream;
    _numMassIRS & stream;
    _numMassVA & stream;
    _numClinical_Microscopy & stream;
    _numAddedToCohort & stream;
    _numRemovedFromCohort & stream;
    _numMDAs & stream;
    _numNmfDeaths & stream;
    _numAntibioticTreatments & stream;
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
  
  // data, per AgeGroup:
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
  vector<double> _inoculationsPerAgeGroup;
  vector<int> _numClinical_FirstDayDeaths;
  vector<int> _numClinical_HospitalFirstDayDeaths;
  vector<int> _numNewInfections;
  vector<int> _numMassITNs;
  vector<int> _numEPI_ITNs;
  vector<int> _numMassIRS;
  vector<int> _numMassVA;
  vector<int> _numAddedToCohort;
  vector<int> _numRemovedFromCohort;
  vector<int> _numMDAs;
  vector<int> _numNmfDeaths;
  vector<int> _numAntibioticTreatments;
  
    // data, per vector species:
    map<string,double> data_Vector_Nv0;
    map<string,double> data_Vector_Nv;
    map<string,double> data_Vector_Ov;
    map<string,double> data_Vector_Sv;
    double data_Vector_EIR_Input;
    double data_Vector_EIR_Simulated;
    
    int _numClinical_RDTs;
    map<string,double> _sumClinical_DrugUsage;
    map<string,double> _sumClinical_DrugUsageIV;
    int _numClinical_Microscopy;
    
  friend class SurveysType;
};

/** Line end character. Use Unix line endings to save a little size. */
const char lineEnd = '\n';

} }
#endif
