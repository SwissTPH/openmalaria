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
    
// -----  Static members  -----

double SurveyAgeGroup::_lowerbound;
vector<double> SurveyAgeGroup::_upperbound;
bitset<NUM_SURVEY_OPTIONS> Survey::active;
bool Survey::_assimilatorMode;


class SurveyCodeMap {
    // Lookup table to translate the strings used in the XML file to the internal enumerated values:
    map<string,SurveyCodes> codeMap;
    
    public:
	SurveyCodeMap () {
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
	    codeMap["innoculationsPerDayOfYear"] = innoculationsPerDayOfYear;
	    codeMap["kappaPerDayOfYear"] = kappaPerDayOfYear;
	    codeMap["innoculationsPerAgeGroup"] = innoculationsPerAgeGroup;
	}
	
	SurveyCodes operator[] (const string s) {
	    map<string,SurveyCodes>::iterator codeIt = codeMap.find (s);
	    if (codeIt == codeMap.end()) {
		ostringstream msg;
		msg << "Unrecognised option: ";
		msg << s;
		throw util::xml_scenario_error(msg.str());
	    }
	    return codeIt->second;
	}
	// reverse-lookup in map; only used for error/debug printing so efficiency is unimportant
	// doesn't ensure code is unique in the map either
	string toString (const SurveyCodes code) {
	    for (map<string,SurveyCodes>::iterator codeIt = codeMap.begin(); codeIt != codeMap.end(); ++codeIt) {
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
    SurveyCodeMap codeMap;
    
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

    for (size_t i = 0;i < numAgeGroups - 1; i++) {
	_upperbound[i] = groups[i].getUpperbound();
    }
    _upperbound[numAgeGroups-1] = DBL_MAX;
}

SurveyAgeGroup::SurveyAgeGroup (double ageYears) {
  if (ageYears < _lowerbound)
    _i = _upperbound.size() - 1;
  else {
    _i = 0;
    while (ageYears > _upperbound[_i])
      _i++;
  }
}


// -----  Non-static members  -----

void Survey::reportTreatment (SurveyAgeGroup ageGroup, int regimen)
{
  switch (regimen) {
    case 1:
      _numTreatments1[ageGroup.i()]++;
      break;
    case 2:
      _numTreatments2[ageGroup.i()]++;
      break;
    case 3:
      _numTreatments3[ageGroup.i()]++;
      break;
    default:
      throw invalid_argument ("Unsupported regimen");
  }
}


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
}


void Survey::writeSummaryArrays (ostream& outputFile, int survey)
{
  if (active[nHost]) {
    writeArray (outputFile, nHost, _assimilatorMode, survey, _numHosts);
  }
  if (active[nInfect]) {
    writeArray (outputFile, nInfect, _assimilatorMode, survey, _numInfectedHosts);
  }
  if (active[nExpectd]) {
    writeArray (outputFile, nExpectd, _assimilatorMode, survey, _numExpectedInfected);
  }
  if (active[nPatent]) {
    writeArray (outputFile, nPatent, _assimilatorMode, survey, _numPatentHosts);
  }
  if (active[sumLogPyrogenThres]) {
    writeArray (outputFile, sumLogPyrogenThres, _assimilatorMode, survey, _sumLogPyrogenicThreshold);
  }
  if (active[sumlogDens]) {
    writeArray (outputFile, sumlogDens, _assimilatorMode, survey, _sumLogDensity);
  }
  if (active[totalInfs]) {
    writeArray (outputFile, totalInfs, _assimilatorMode, survey, _sumInfections);
  }
  if (active[nTransmit]) {
    writeArray (outputFile, nTransmit, _assimilatorMode, survey, _numTransmittingHosts);
  }
  if (active[totalPatentInf]) {
    writeArray (outputFile, totalPatentInf, _assimilatorMode, survey, _sumPatentInfections);
  }
  if (active[sumPyrogenThresh]) {
    writeArray (outputFile, sumPyrogenThresh, _assimilatorMode, survey, _sumPyrogenicThreshold);
  }
  if (active[nTreatments1]) {
    writeArray (outputFile, nTreatments1, _assimilatorMode, survey, _numTreatments1);
  }
  if (active[nTreatments2]) {
    writeArray (outputFile, nTreatments2, _assimilatorMode, survey, _numTreatments2);
  }
  if (active[nTreatments3]) {
    writeArray (outputFile, nTreatments3, _assimilatorMode, survey, _numTreatments3);
  }
  if (active[nUncomp]) {
    writeArray (outputFile, nUncomp, _assimilatorMode, survey, _numUncomplicatedEpisodes);
  }
  if (active[nSevere]) {
    writeArray (outputFile, nSevere, _assimilatorMode, survey, _numSevereEpisodes);
  }
  if (active[nSeq]) {
    writeArray (outputFile, nSeq, _assimilatorMode, survey, _numSequelae);
  }
  if (active[nHospitalDeaths]) {
    writeArray (outputFile, nHospitalDeaths, _assimilatorMode, survey, _numHospitalDeaths);
  }
  if (active[nIndDeaths]) {
    writeArray (outputFile, nIndDeaths, _assimilatorMode, survey, _numIndirectDeaths);
  }
  if (active[nDirDeaths]) {
    writeArray (outputFile, nDirDeaths, _assimilatorMode, survey, _numDirectDeaths);
  }
  if (active[nEPIVaccinations]) {
    writeArray (outputFile, nEPIVaccinations, _assimilatorMode, survey, _numEPIVaccinations);
  }
  if (active[nMassVaccinations]) {
    writeArray (outputFile, nMassVaccinations, _assimilatorMode, survey, _numMassVaccinations);
  }
  if (active[nHospitalRecovs]) {
    writeArray (outputFile, nHospitalRecovs, _assimilatorMode, survey, _numHospitalRecoveries);
  }
  if (active[nHospitalSeqs]) {
    writeArray (outputFile, nHospitalSeqs, _assimilatorMode, survey, _numHospitalSequelae);
  }
  if (active[nIPTDoses]) {
    writeArray (outputFile, nIPTDoses, _assimilatorMode, survey, _numIPTDoses);
  }
  if (active[annAvgK]) {
    writeArray (outputFile, annAvgK, _assimilatorMode, survey, _annualAverageKappa);
  }
  if (active[nNMFever]) {
    writeArray (outputFile, nNMFever, _assimilatorMode, survey, _numNonMalariaFevers);
  }

  if (active[innoculationsPerDayOfYear]) {
    writeArray (outputFile, innoculationsPerDayOfYear, _assimilatorMode, survey, _innoculationsPerDayOfYear);
  }
  if (active[kappaPerDayOfYear]) {
    writeArray (outputFile, kappaPerDayOfYear, _assimilatorMode, survey, _kappaPerDayOfYear);
  }
  if (active[innoculationsPerAgeGroup]) {
    writeArray (outputFile, innoculationsPerAgeGroup, _assimilatorMode, survey, _innoculationsPerAgeGroup);
  }
}


template <class T>
void writeArray (ostream& file, int measure, bool assimilatorMode, int survey, vector<T>& array)
{
  for (int j = 0; j < (int) array.size() - 1; j++) { // Don't write out last age-group
    if (!assimilatorMode)
      file << survey << "\t" << j + 1 << "\t" << measure;
    file << "\t" << array[j] << lineEnd;
  }
}

template <class T>
void writeArray (ostream& file, int measure, bool assimilatorMode, int survey, T& value)
{
  if (!assimilatorMode)
    file << survey << "\t" << 0 << "\t" << measure;
  file << "\t" << value << lineEnd;
}

}