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

#include "Monitoring/Survey.h"
#include "inputData.h"
#include "util/errors.hpp"

#include <stdexcept>

namespace OM { namespace Monitoring {
 
// -----  Utility forward declarations  -----

template <class T>
void writeValue(ostream& file, int measure, bool assimilationMode, int survey, T& value);

void writeMap (ostream& file, int measure, bool assimilatorMode, int survey, map<string,double>& data);

template <class T>
void writePerAgeGroup(ostream& file, int measure, bool assimilationMode, int survey, vector<T>& array);


// -----  Static members  -----

double AgeGroup::_lowerbound;
vector<double> AgeGroup::_upperbound;
bitset<SM::NUM_SURVEY_OPTIONS> Survey::active;
bool Survey::_assimilatorMode;


class SurveyMeasureMap {
    // Lookup table to translate the strings used in the XML file to the internal enumerated values:
    map<string,SM::SurveyMeasure> codeMap;
    
    public:
	SurveyMeasureMap () {
	    codeMap["nHost"] = SM::nHost;
	    codeMap["nInfect"] = SM::nInfect;
	    codeMap["nExpectd"] = SM::nExpectd;
	    codeMap["nPatent"] = SM::nPatent;
	    codeMap["sumLogPyrogenThres"] = SM::sumLogPyrogenThres;
	    codeMap["sumlogDens"] = SM::sumlogDens;
	    codeMap["totalInfs"] = SM::totalInfs;
	    codeMap["nTransmit"] = SM::nTransmit;
	    codeMap["totalPatentInf"] = SM::totalPatentInf;
	    codeMap["contrib"] = SM::contrib;
	    codeMap["sumPyrogenThresh"] = SM::sumPyrogenThresh;
	    codeMap["nTreatments1"] = SM::nTreatments1;
	    codeMap["nTreatments2"] = SM::nTreatments2;
	    codeMap["nTreatments3"] = SM::nTreatments3;
	    codeMap["nUncomp"] = SM::nUncomp;
	    codeMap["nSevere"] = SM::nSevere;
	    codeMap["nSeq"] = SM::nSeq;
	    codeMap["nHospitalDeaths"] = SM::nHospitalDeaths;
	    codeMap["nIndDeaths"] = SM::nIndDeaths;
	    codeMap["nDirDeaths"] = SM::nDirDeaths;
	    codeMap["nEPIVaccinations"] = SM::nEPIVaccinations;
	    codeMap["allCauseIMR"] = SM::allCauseIMR;
	    codeMap["nMassVaccinations"] = SM::nMassVaccinations;
	    codeMap["nHospitalRecovs"] = SM::nHospitalRecovs;
	    codeMap["nHospitalSeqs"] = SM::nHospitalSeqs;
	    codeMap["nIPTDoses"] = SM::nIPTDoses;
	    codeMap["annAvgK"] = SM::annAvgK;
	    codeMap["nNMFever"] = SM::nNMFever;
	    codeMap["innoculationsPerAgeGroup"] = SM::innoculationsPerAgeGroup;
	    codeMap["Vector_Nv0"] = SM::Vector_Nv0;
	    codeMap["Vector_Nv"] = SM::Vector_Nv;
	    codeMap["Vector_Ov"] = SM::Vector_Ov;
	    codeMap["Vector_Sv"] = SM::Vector_Sv;
	    codeMap["Vector_EIR_Input"] = SM::Vector_EIR_Input;
	    codeMap["Vector_EIR_Simulated"] = SM::Vector_EIR_Simulated;
	    codeMap["Clinical_RDTs"] = SM::Clinical_RDTs;
	    codeMap["Clinical_DrugUsage"] = SM::Clinical_DrugUsage;
	    codeMap["Clinical_FirstDayDeaths"] = SM::Clinical_FirstDayDeaths;
	    codeMap["Clinical_HospitalFirstDayDeaths"] = SM::Clinical_HospitalFirstDayDeaths;
	    codeMap["nNewInfections"] = SM::nNewInfections;
	    codeMap["nMassITNs"] = SM::nMassITNs;
	    codeMap["nEPI_ITNs"] = SM::nEPI_ITNs;
	    codeMap["nMassIRS"] = SM::nMassIRS;
	    codeMap["nMassVA"] = SM::nMassVA;
	    codeMap["Clinical_Microscopy"] = SM::Clinical_Microscopy;
	}
	
	SM::SurveyMeasure operator[] (const string s) {
	    map<string,SM::SurveyMeasure>::iterator codeIt = codeMap.find (s);
	    if (codeIt == codeMap.end()) {
		ostringstream msg;
		msg << "Unrecognised survey option: \"";
		msg << s << '"';
		throw util::xml_scenario_error(msg.str());
	    }
	    return codeIt->second;
	}
	// reverse-lookup in map; only used for error/debug printing so efficiency is unimportant
	// doesn't ensure code is unique in the map either
	string toString (const SM::SurveyMeasure code) {
	    for (map<string,SM::SurveyMeasure>::iterator codeIt = codeMap.begin(); codeIt != codeMap.end(); ++codeIt) {
		if (codeIt->second == code)
		    return codeIt->first;
	    }
	    throw runtime_error ("toString called with unknown code");	// this is a code error
	}
};


void Survey::init () {
    AgeGroup::init ();
    
    // by default, none are active
    active.reset ();
    SurveyMeasureMap codeMap;
    
    scnXml::OptionSet::OptionSequence sOSeq = InputData().getMonitoring().getSurveyOptions().getOption();
    for (scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
	active[codeMap[it->getName()]] = it->getValue();
    }
    
    _assimilatorMode = InputData().getAssimMode();
}
void AgeGroup::init () {
    const scnXml::Monitoring& mon = InputData().getMonitoring();
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
    _upperbound[numAgeGroups-1] = numeric_limits<double>::infinity();
}

void AgeGroup::update (double ageYears) {
    while (ageYears > _upperbound[_i])
	_i++;
}


// -----  Non-static members  -----


void Survey::allocate ()
{
  size_t numAgeGroups = AgeGroup::getNumGroups();
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
    
    _numClinical_RDTs = 0;
    _numClinical_FirstDayDeaths.resize (numAgeGroups);
    _numClinical_HospitalFirstDayDeaths.resize (numAgeGroups);
    _numNewInfections.resize (numAgeGroups);
    _numMassITNs.resize (numAgeGroups);
    _numEPI_ITNs.resize (numAgeGroups);
    _numMassIRS.resize (numAgeGroups);
    _numMassVA.resize (numAgeGroups);
}


void Survey::writeSummaryArrays (ostream& outputFile, int survey)
{
  if (active[SM::nHost]) {
    writePerAgeGroup (outputFile, SM::nHost, _assimilatorMode, survey, _numHosts);
  }
  if (active[SM::nInfect]) {
    writePerAgeGroup (outputFile, SM::nInfect, _assimilatorMode, survey, _numInfectedHosts);
  }
  if (active[SM::nExpectd]) {
    writePerAgeGroup (outputFile, SM::nExpectd, _assimilatorMode, survey, _numExpectedInfected);
  }
  if (active[SM::nPatent]) {
    writePerAgeGroup (outputFile, SM::nPatent, _assimilatorMode, survey, _numPatentHosts);
  }
  if (active[SM::sumLogPyrogenThres]) {
    writePerAgeGroup (outputFile, SM::sumLogPyrogenThres, _assimilatorMode, survey, _sumLogPyrogenicThreshold);
  }
  if (active[SM::sumlogDens]) {
    writePerAgeGroup (outputFile, SM::sumlogDens, _assimilatorMode, survey, _sumLogDensity);
  }
  if (active[SM::totalInfs]) {
    writePerAgeGroup (outputFile, SM::totalInfs, _assimilatorMode, survey, _sumInfections);
  }
  if (active[SM::nTransmit]) {
    writeValue (outputFile, SM::nTransmit, _assimilatorMode, survey, _numTransmittingHosts);
  }
  if (active[SM::totalPatentInf]) {
    writePerAgeGroup (outputFile, SM::totalPatentInf, _assimilatorMode, survey, _sumPatentInfections);
  }
  if (active[SM::sumPyrogenThresh]) {
    writePerAgeGroup (outputFile, SM::sumPyrogenThresh, _assimilatorMode, survey, _sumPyrogenicThreshold);
  }
  if (active[SM::nTreatments1]) {
    writePerAgeGroup (outputFile, SM::nTreatments1, _assimilatorMode, survey, _numTreatments1);
  }
  if (active[SM::nTreatments2]) {
    writePerAgeGroup (outputFile, SM::nTreatments2, _assimilatorMode, survey, _numTreatments2);
  }
  if (active[SM::nTreatments3]) {
    writePerAgeGroup (outputFile, SM::nTreatments3, _assimilatorMode, survey, _numTreatments3);
  }
  if (active[SM::nUncomp]) {
    writePerAgeGroup (outputFile, SM::nUncomp, _assimilatorMode, survey, _numUncomplicatedEpisodes);
  }
  if (active[SM::nSevere]) {
    writePerAgeGroup (outputFile, SM::nSevere, _assimilatorMode, survey, _numSevereEpisodes);
  }
  if (active[SM::nSeq]) {
    writePerAgeGroup (outputFile, SM::nSeq, _assimilatorMode, survey, _numSequelae);
  }
  if (active[SM::nHospitalDeaths]) {
    writePerAgeGroup (outputFile, SM::nHospitalDeaths, _assimilatorMode, survey, _numHospitalDeaths);
  }
  if (active[SM::nIndDeaths]) {
    writePerAgeGroup (outputFile, SM::nIndDeaths, _assimilatorMode, survey, _numIndirectDeaths);
  }
  if (active[SM::nDirDeaths]) {
    writePerAgeGroup (outputFile, SM::nDirDeaths, _assimilatorMode, survey, _numDirectDeaths);
  }
  if (active[SM::nEPIVaccinations]) {
    writePerAgeGroup (outputFile, SM::nEPIVaccinations, _assimilatorMode, survey, _numEPIVaccinations);
  }
  if (active[SM::nMassVaccinations]) {
    writePerAgeGroup (outputFile, SM::nMassVaccinations, _assimilatorMode, survey, _numMassVaccinations);
  }
  if (active[SM::nHospitalRecovs]) {
    writePerAgeGroup (outputFile, SM::nHospitalRecovs, _assimilatorMode, survey, _numHospitalRecoveries);
  }
  if (active[SM::nHospitalSeqs]) {
    writePerAgeGroup (outputFile, SM::nHospitalSeqs, _assimilatorMode, survey, _numHospitalSequelae);
  }
  if (active[SM::nIPTDoses]) {
    writePerAgeGroup (outputFile, SM::nIPTDoses, _assimilatorMode, survey, _numIPTDoses);
  }
  if (active[SM::annAvgK]) {
    writeValue (outputFile, SM::annAvgK, _assimilatorMode, survey, _annualAverageKappa);
  }
  if (active[SM::nNMFever]) {
    writePerAgeGroup (outputFile, SM::nNMFever, _assimilatorMode, survey, _numNonMalariaFevers);
  }

  if (active[SM::innoculationsPerAgeGroup]) {
    writePerAgeGroup (outputFile, SM::innoculationsPerAgeGroup, _assimilatorMode, survey, _innoculationsPerAgeGroup);
  }
  
  if (active[SM::Vector_Nv0]) {
    writeMap (outputFile, SM::Vector_Nv0, _assimilatorMode, survey, data_Vector_Nv0);
  }
  if (active[SM::Vector_Nv]) {
    writeMap (outputFile, SM::Vector_Nv, _assimilatorMode, survey, data_Vector_Nv);
  }
  if (active[SM::Vector_Ov]) {
    writeMap (outputFile, SM::Vector_Ov, _assimilatorMode, survey, data_Vector_Ov);
  }
  if (active[SM::Vector_Sv]) {
    writeMap (outputFile, SM::Vector_Sv, _assimilatorMode, survey, data_Vector_Sv);
  }
  if (active[SM::Vector_EIR_Input]) {
    writeValue (outputFile, SM::Vector_EIR_Input, _assimilatorMode, survey, data_Vector_EIR_Input);
  }
  if (active[SM::Vector_EIR_Simulated]) {
    writeValue (outputFile, SM::Vector_EIR_Simulated, _assimilatorMode, survey, data_Vector_EIR_Simulated);
  }
  if (active[SM::Clinical_RDTs]) {
      writeValue (outputFile, SM::Clinical_RDTs, _assimilatorMode, survey, _numClinical_RDTs);
  }
  if (active[SM::Clinical_DrugUsage]) {
      writeMap (outputFile, SM::Clinical_DrugUsage, _assimilatorMode, survey, _sumClinical_DrugUsage);
  }
  if (active[SM::Clinical_FirstDayDeaths]) {
      writePerAgeGroup (outputFile, SM::Clinical_FirstDayDeaths, _assimilatorMode, survey, _numClinical_FirstDayDeaths);
  }
  if (active[SM::Clinical_HospitalFirstDayDeaths]) {
      writePerAgeGroup (outputFile, SM::Clinical_HospitalFirstDayDeaths, _assimilatorMode, survey, _numClinical_HospitalFirstDayDeaths);
  }
  if (active[SM::nNewInfections]) {
      writePerAgeGroup (outputFile, SM::nNewInfections, _assimilatorMode, survey, _numNewInfections);
  }
  if (active[SM::nMassITNs]) {
      writePerAgeGroup (outputFile, SM::nMassITNs, _assimilatorMode, survey, _numMassITNs);
  }
  if (active[SM::nEPI_ITNs]) {
      writePerAgeGroup (outputFile, SM::nEPI_ITNs, _assimilatorMode, survey, _numEPI_ITNs);
  }
  if (active[SM::nMassIRS]) {
      writePerAgeGroup (outputFile, SM::nMassIRS, _assimilatorMode, survey, _numMassIRS);
  }
  if (active[SM::nMassVA]) {
      writePerAgeGroup (outputFile, SM::nMassVA, _assimilatorMode, survey, _numMassVA);
  }
  if (active[SM::Clinical_Microscopy]) {
      writeValue (outputFile, SM::Clinical_Microscopy, _assimilatorMode, survey, _numClinical_Microscopy);
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

} }