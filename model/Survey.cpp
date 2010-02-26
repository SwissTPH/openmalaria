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

#include "Survey.h"
#include "inputData.h"
#include "util/errors.hpp"

#include <cfloat>
#include <stdexcept>

namespace OM {
 
// -----  Utility forward declarations  -----

template <class T>
void writeValue(ostream& file, int measure, bool assimilationMode, int survey, T& value);

void writeMap (ostream& file, int measure, bool assimilatorMode, int survey, map<string,double>& data);

template <class T>
void writePerAgeGroup(ostream& file, int measure, bool assimilationMode, int survey, vector<T>& array);


// -----  Static members  -----

double SurveyAgeGroup::_lowerbound;
vector<double> SurveyAgeGroup::_upperbound;
bitset<NUM_SURVEY_OPTIONS> Survey::active;
bool Survey::_assimilatorMode;


class SurveyMeasureMap {
    // Lookup table to translate the strings used in the XML file to the internal enumerated values:
    map<string,SurveyMeasure> codeMap;
    
    public:
	SurveyMeasureMap () {
	    codeMap["nHost"] = nHost;
	    codeMap["nInfect"] = nInfect;
	    codeMap["nExpectd"] = nExpectd;
	    codeMap["nPatent"] = nPatent;
	    codeMap["sumLogPyrogenThres"] = sumLogPyrogenThres;
	    codeMap["sumlogDens"] = sumlogDens;
	    codeMap["totalInfs"] = totalInfs;
	    codeMap["nTransmit"] = nTransmit;
	    codeMap["totalPatentInf"] = totalPatentInf;
	    codeMap["contrib"] = contrib;
	    codeMap["sumPyrogenThresh"] = sumPyrogenThresh;
	    codeMap["nTreatments1"] = nTreatments1;
	    codeMap["nTreatments2"] = nTreatments2;
	    codeMap["nTreatments3"] = nTreatments3;
	    codeMap["nUncomp"] = nUncomp;
	    codeMap["nSevere"] = nSevere;
	    codeMap["nSeq"] = nSeq;
	    codeMap["nHospitalDeaths"] = nHospitalDeaths;
	    codeMap["nIndDeaths"] = nIndDeaths;
	    codeMap["nDirDeaths"] = nDirDeaths;
	    codeMap["nEPIVaccinations"] = nEPIVaccinations;
	    codeMap["imr_summary"] = imr_summary;
	    codeMap["nMassVaccinations"] = nMassVaccinations;
	    codeMap["nHospitalRecovs"] = nHospitalRecovs;
	    codeMap["nHospitalSeqs"] = nHospitalSeqs;
	    codeMap["nIPTDoses"] = nIPTDoses;
	    codeMap["annAvgK"] = annAvgK;
	    codeMap["nNMFever"] = nNMFever;
	    codeMap["innoculationsPerAgeGroup"] = innoculationsPerAgeGroup;
	    codeMap["Vector_Nv0"] = Vector_Nv0;
	    codeMap["Vector_Nv"] = Vector_Nv;
	    codeMap["Vector_Ov"] = Vector_Ov;
	    codeMap["Vector_Sv"] = Vector_Sv;
	    codeMap["Vector_EIR_Input"] = Vector_EIR_Input;
	    codeMap["Vector_EIR_Simulated"] = Vector_EIR_Simulated;
	}
	
	SurveyMeasure operator[] (const string s) {
	    map<string,SurveyMeasure>::iterator codeIt = codeMap.find (s);
	    if (codeIt == codeMap.end()) {
		ostringstream msg;
		msg << "Unrecognised survey option: ";
		msg << s;
		throw util::xml_scenario_error(msg.str());
	    }
	    return codeIt->second;
	}
	// reverse-lookup in map; only used for error/debug printing so efficiency is unimportant
	// doesn't ensure code is unique in the map either
	string toString (const SurveyMeasure code) {
	    for (map<string,SurveyMeasure>::iterator codeIt = codeMap.begin(); codeIt != codeMap.end(); ++codeIt) {
		if (codeIt->second == code)
		    return codeIt->first;
	    }
	    throw runtime_error ("toString called with unknown code");	// this is a code error
	}
};


void Survey::init () {
    SurveyAgeGroup::init ();
    
    // by default, none are active
    active.reset ();
    SurveyMeasureMap codeMap;
    
    scnXml::OptionSet::OptionSequence sOSeq = InputData.getMonitoring().getSurveyOptions().getOption();
    for (scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
	active[codeMap[it->getName()]] = it->getValue();
    }
    
    _assimilatorMode = InputData.getScenario().getAssimMode();
}
void SurveyAgeGroup::init () {
    const scnXml::Monitoring& mon = InputData.getMonitoring();
    const scnXml::AgeGroup::GroupSequence& groups = mon.getAgeGroup().getGroup();
    /* note that the last age group includes individuals who are        *
    * either younger than Lowerbound or older than the last Upperbound */
    size_t numAgeGroups = groups.size() + 1;
    _upperbound.resize (numAgeGroups);
    _lowerbound = mon.getAgeGroup().getLowerbound();
    if (!(_lowerbound <= 0.0))
	throw util::xml_scenario_error ("Expected survey age-group lowerbound of 0");

    for (size_t i = 0;i < numAgeGroups - 1; i++) {
	_upperbound[i] = groups[i].getUpperbound();
    }
    _upperbound[numAgeGroups-1] = DBL_MAX;
}

void SurveyAgeGroup::update (double ageYears) {
    while (ageYears > _upperbound[_i])
	_i++;
}


// -----  Non-static members  -----


void Survey::allocate ()
{
  size_t numAgeGroups = SurveyAgeGroup::getNumGroups();
  _numHosts.resize (numAgeGroups);
  _numInfectedHosts.resize (numAgeGroups);
  _numExpectedInfected.resize (numAgeGroups);
  _numPatentHosts.resize (numAgeGroups);
  _sumLogPyrogenicThreshold.resize (numAgeGroups);
  _sumLogDensity.resize (numAgeGroups);
  _sumInfections.resize (numAgeGroups);
  _numTransmittingHosts = numeric_limits<double>::signaling_NaN();
  _sumPatentInfections.resize (numAgeGroups);
  _sumPyrogenicThreshold.resize (numAgeGroups);
  _numTreatments1.resize (numAgeGroups);
  _numTreatments2.resize (numAgeGroups);
  _numTreatments3.resize (numAgeGroups);
  _numUncomplicatedEpisodes.resize (numAgeGroups);
  _numSevereEpisodes.resize (numAgeGroups);
  _numSequelae.resize (numAgeGroups);
  _numHospitalDeaths.resize (numAgeGroups);
  _numIndirectDeaths.resize (numAgeGroups);
  _numDirectDeaths.resize (numAgeGroups);
  _numEPIVaccinations.resize (numAgeGroups);
  _numMassVaccinations.resize (numAgeGroups);
  _numHospitalRecoveries.resize (numAgeGroups);
  _numHospitalSequelae.resize (numAgeGroups);
  _numIPTDoses.resize (numAgeGroups);
  _annualAverageKappa = numeric_limits<double>::signaling_NaN();
  _numNonMalariaFevers.resize (numAgeGroups);
  // These 3 are set as a whole and so don't require allocation:
//   _innoculationsPerDayOfYear;
//   _kappaPerDayOfYear;
//   _innoculationsPerAgeGroup;
    
    data_Vector_EIR_Input =numeric_limits<double>::signaling_NaN() ;
    data_Vector_EIR_Simulated = numeric_limits<double>::signaling_NaN();
    
    data_Clinical_HospitalEntries = 0;
    data_Clinical_HospitalizationDays = 0;
    data_Clinical_RDTs = 0;
}


void Survey::writeSummaryArrays (ostream& outputFile, int survey)
{
  if (active[nHost]) {
    writePerAgeGroup (outputFile, nHost, _assimilatorMode, survey, _numHosts);
  }
  if (active[nInfect]) {
    writePerAgeGroup (outputFile, nInfect, _assimilatorMode, survey, _numInfectedHosts);
  }
  if (active[nExpectd]) {
    writePerAgeGroup (outputFile, nExpectd, _assimilatorMode, survey, _numExpectedInfected);
  }
  if (active[nPatent]) {
    writePerAgeGroup (outputFile, nPatent, _assimilatorMode, survey, _numPatentHosts);
  }
  if (active[sumLogPyrogenThres]) {
    writePerAgeGroup (outputFile, sumLogPyrogenThres, _assimilatorMode, survey, _sumLogPyrogenicThreshold);
  }
  if (active[sumlogDens]) {
    writePerAgeGroup (outputFile, sumlogDens, _assimilatorMode, survey, _sumLogDensity);
  }
  if (active[totalInfs]) {
    writePerAgeGroup (outputFile, totalInfs, _assimilatorMode, survey, _sumInfections);
  }
  if (active[nTransmit]) {
    writeValue (outputFile, nTransmit, _assimilatorMode, survey, _numTransmittingHosts);
  }
  if (active[totalPatentInf]) {
    writePerAgeGroup (outputFile, totalPatentInf, _assimilatorMode, survey, _sumPatentInfections);
  }
  if (active[sumPyrogenThresh]) {
    writePerAgeGroup (outputFile, sumPyrogenThresh, _assimilatorMode, survey, _sumPyrogenicThreshold);
  }
  if (active[nTreatments1]) {
    writePerAgeGroup (outputFile, nTreatments1, _assimilatorMode, survey, _numTreatments1);
  }
  if (active[nTreatments2]) {
    writePerAgeGroup (outputFile, nTreatments2, _assimilatorMode, survey, _numTreatments2);
  }
  if (active[nTreatments3]) {
    writePerAgeGroup (outputFile, nTreatments3, _assimilatorMode, survey, _numTreatments3);
  }
  if (active[nUncomp]) {
    writePerAgeGroup (outputFile, nUncomp, _assimilatorMode, survey, _numUncomplicatedEpisodes);
  }
  if (active[nSevere]) {
    writePerAgeGroup (outputFile, nSevere, _assimilatorMode, survey, _numSevereEpisodes);
  }
  if (active[nSeq]) {
    writePerAgeGroup (outputFile, nSeq, _assimilatorMode, survey, _numSequelae);
  }
  if (active[nHospitalDeaths]) {
    writePerAgeGroup (outputFile, nHospitalDeaths, _assimilatorMode, survey, _numHospitalDeaths);
  }
  if (active[nIndDeaths]) {
    writePerAgeGroup (outputFile, nIndDeaths, _assimilatorMode, survey, _numIndirectDeaths);
  }
  if (active[nDirDeaths]) {
    writePerAgeGroup (outputFile, nDirDeaths, _assimilatorMode, survey, _numDirectDeaths);
  }
  if (active[nEPIVaccinations]) {
    writePerAgeGroup (outputFile, nEPIVaccinations, _assimilatorMode, survey, _numEPIVaccinations);
  }
  if (active[nMassVaccinations]) {
    writePerAgeGroup (outputFile, nMassVaccinations, _assimilatorMode, survey, _numMassVaccinations);
  }
  if (active[nHospitalRecovs]) {
    writePerAgeGroup (outputFile, nHospitalRecovs, _assimilatorMode, survey, _numHospitalRecoveries);
  }
  if (active[nHospitalSeqs]) {
    writePerAgeGroup (outputFile, nHospitalSeqs, _assimilatorMode, survey, _numHospitalSequelae);
  }
  if (active[nIPTDoses]) {
    writePerAgeGroup (outputFile, nIPTDoses, _assimilatorMode, survey, _numIPTDoses);
  }
  if (active[annAvgK]) {
    writeValue (outputFile, annAvgK, _assimilatorMode, survey, _annualAverageKappa);
  }
  if (active[nNMFever]) {
    writePerAgeGroup (outputFile, nNMFever, _assimilatorMode, survey, _numNonMalariaFevers);
  }

  if (active[innoculationsPerAgeGroup]) {
    writePerAgeGroup (outputFile, innoculationsPerAgeGroup, _assimilatorMode, survey, _innoculationsPerAgeGroup);
  }
  
  if (active[Vector_Nv0]) {
    writeMap (outputFile, Vector_Nv0, _assimilatorMode, survey, data_Vector_Nv0);
  }
  if (active[Vector_Nv]) {
    writeMap (outputFile, Vector_Nv, _assimilatorMode, survey, data_Vector_Nv);
  }
  if (active[Vector_Ov]) {
    writeMap (outputFile, Vector_Ov, _assimilatorMode, survey, data_Vector_Ov);
  }
  if (active[Vector_Sv]) {
    writeMap (outputFile, Vector_Sv, _assimilatorMode, survey, data_Vector_Sv);
  }
  if (active[Vector_EIR_Input]) {
    writeValue (outputFile, Vector_EIR_Input, _assimilatorMode, survey, data_Vector_EIR_Input);
  }
  if (active[Vector_EIR_Simulated]) {
    writeValue (outputFile, Vector_EIR_Simulated, _assimilatorMode, survey, data_Vector_EIR_Simulated);
  }
  if (active[Clinical_HospitalEntries]) {
      writeValue (outputFile, Clinical_HospitalEntries, _assimilatorMode, survey, data_Clinical_HospitalEntries);
  }
  if (active[Clinical_HospitalizationDays]) {
      writeValue (outputFile, Clinical_HospitalizationDays, _assimilatorMode, survey, data_Clinical_HospitalizationDays);
  }
  if (active[Clinical_RDTs]) {
      writeValue (outputFile, Clinical_RDTs, _assimilatorMode, survey, data_Clinical_RDTs);
  }
  if (active[Clinical_DrugUsage]) {
      writeMap (outputFile, Clinical_DrugUsage, _assimilatorMode, survey, data_Clinical_DrugUsage);
  }
}


template <class T>
void writeValue (ostream& file, int measure, bool assimilatorMode, int survey, T& value)
{
  if (!assimilatorMode)
    file << survey << "\t" << 0 << "\t" << measure;
  file << "\t" << value << lineEnd;
}

void writeMap (ostream& file, int measure, bool assimilatorMode, int survey, map<string,double>& data)
{
    for (map<string,double>::const_iterator it = data.begin(); it != data.end(); ++it) {
	if (!assimilatorMode)
	    file << survey << "\t" << it->first << "\t" << measure;
	file << "\t" << it->second << lineEnd;
    }
}

template <class T>
void writePerAgeGroup (ostream& file, int measure, bool assimilatorMode, int survey, vector<T>& array)
{
  for (int j = 0; j < (int) array.size() - 1; j++) { // Don't write out last age-group
    if (!assimilatorMode)
      file << survey << "\t" << j + 1 << "\t" << measure;
    file << "\t" << array[j] << lineEnd;
  }
}

}