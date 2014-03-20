/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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
#include "util/errors.h"
#include "Host/Human.h"
#include "schema/monitoring.h"

#include <boost/static_assert.hpp>
#include <stdexcept>

namespace OM { namespace Monitoring {
 
// -----  Utility forward declarations  -----

template <class T>
void writeValue(ostream& file, int measure, int survey, T& value);

void writeMap (ostream& file, int measure, int survey, map<string,double>& data);

size_t intReportMappings[Report::MI_NUM];
size_t dblReportMappings[Report::MD_NUM];

// -----  Static members  -----

double AgeGroup::_lowerbound;
vector<double> AgeGroup::_upperbound;
bitset<SM::NUM_SURVEY_OPTIONS> Survey::active;


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
        codeMap["inputEIR"] = SM::inputEIR;
        codeMap["simulatedEIR"] = SM::simulatedEIR;
        codeMap["Clinical_RDTs"] = SM::Clinical_RDTs;
        codeMap["Clinical_DrugUsage"] = SM::Clinical_DrugUsage;
        codeMap["Clinical_DrugUsageIV"] = SM::Clinical_DrugUsageIV;
        codeMap["Clinical_FirstDayDeaths"] = SM::Clinical_FirstDayDeaths;
        codeMap["Clinical_HospitalFirstDayDeaths"] = SM::Clinical_HospitalFirstDayDeaths;
        codeMap["nNewInfections"] = SM::nNewInfections;
        codeMap["nMassITNs"] = SM::nMassITNs;
        codeMap["nEPI_ITNs"] = SM::nEPI_ITNs;
        codeMap["nMassIRS"] = SM::nMassIRS;
        codeMap["nMassGVI"] = SM::nMassGVI;
        codeMap["Clinical_Microscopy"] = SM::Clinical_Microscopy;
        codeMap["nAddedToCohort"] = SM::nAddedToCohort;
        codeMap["nRemovedFromCohort"] = SM::nRemovedFromCohort;
        codeMap["nMDAs"] = SM::nMDAs;
        codeMap["nMassScreenings"] = SM::nMassScreenings;
        codeMap["nNmfDeaths"] = SM::nNmfDeaths;
        codeMap["nAntibioticTreatments"] = SM::nAntibioticTreatments;
        codeMap["nCtsIRS"] = SM::nCtsIRS;
        codeMap["nCtsGVI"] = SM::nCtsGVI;
        codeMap["nCtsMDA"] = SM::nCtsMDA;
        codeMap["nCtsScreenings"] = SM::nCtsScreenings;
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
        throw TRACED_EXCEPTION_DEFAULT ("toString called with unknown code");	// this is a code error
    }
};


void Survey::init (const scnXml::Monitoring& monitoring) {
    intReportMappings[Report::MI_HOSTS] = SM::nHost;
    intReportMappings[Report::MI_INFECTED_HOSTS] = SM::nInfect;
    intReportMappings[Report::MI_PATENT_HOSTS] = SM::nPatent;
    intReportMappings[Report::MI_INFECTIONS] = SM::totalInfs;
    intReportMappings[Report::MI_PATENT_INFECTIONS] = SM::totalPatentInf;
    intReportMappings[Report::MI_TREATMENTS_1] = SM::nTreatments1;
    intReportMappings[Report::MI_TREATMENTS_2] = SM::nTreatments2;
    intReportMappings[Report::MI_TREATMENTS_3] = SM::nTreatments3;
    intReportMappings[Report::MI_UNCOMPLICATED_EPISODES] = SM::nUncomp;
    intReportMappings[Report::MI_SEVERE_EPISODES] = SM::nSevere;
    intReportMappings[Report::MI_SEQUELAE] = SM::nSeq;
    intReportMappings[Report::MI_HOSPITAL_DEATHS] = SM::nHospitalDeaths;
    intReportMappings[Report::MI_INDIRECT_DEATHS] = SM::nIndDeaths;
    intReportMappings[Report::MI_DIRECT_DEATHS] = SM::nDirDeaths;
    intReportMappings[Report::MI_VACCINATION_TIMED] = SM::nMassVaccinations;
    intReportMappings[Report::MI_VACCINATION_CTS] = SM::nEPIVaccinations;
    intReportMappings[Report::MI_HOSPITAL_RECOVERIES] = SM::nHospitalRecovs;
    intReportMappings[Report::MI_HOSPITAL_SEQUELAE] = SM::nHospitalSeqs;
    intReportMappings[Report::MI_NON_MALARIA_FEVERS] = SM::nNMFever;
    intReportMappings[Report::MI_NEW_INFECTIONS] = SM::nNewInfections;
    intReportMappings[Report::MI_ITN_TIMED] = SM::nMassITNs;
    intReportMappings[Report::MI_ITN_CTS] = SM::nEPI_ITNs;
    intReportMappings[Report::MI_IRS_TIMED] = SM::nMassIRS;
    intReportMappings[Report::MI_IRS_CTS] = SM::nCtsIRS;
    intReportMappings[Report::MI_GVI_TIMED] = SM::nMassGVI;
    intReportMappings[Report::MI_GVI_CTS] = SM::nCtsGVI;
    intReportMappings[Report::MI_MDA_TIMED] = SM::nMDAs;
    intReportMappings[Report::MI_MDA_CTS] = SM::nCtsMDA;
    intReportMappings[Report::MI_SCREENING_TIMED] = SM::nMassScreenings;
    intReportMappings[Report::MI_SCREENING_CTS] = SM::nCtsScreenings;
    intReportMappings[Report::MI_NMF_DEATHS] = SM::nNmfDeaths;
    intReportMappings[Report::MI_NMF_TREATMENTS] = SM::nAntibioticTreatments;
    intReportMappings[Report::MI_FIRST_DAY_DEATHS] = SM::Clinical_FirstDayDeaths;
    intReportMappings[Report::MI_HOSPITAL_FIRST_DAY_DEATHS] = SM::Clinical_HospitalFirstDayDeaths;
    intReportMappings[Report::MI_NUM_ADDED_COHORT] = SM::nAddedToCohort;
    intReportMappings[Report::MI_NUM_REMOVED_COHORT] = SM::nRemovedFromCohort;
    
    dblReportMappings[Report::MD_EXPECTED_INFECTED] = SM::nExpectd;
    dblReportMappings[Report::MD_LOG_PYROGENIC_THRESHOLD] = SM::sumLogPyrogenThres;
    dblReportMappings[Report::MD_LOG_DENSITY] = SM::sumlogDens;
    dblReportMappings[Report::MD_PYROGENIC_THRESHOLD] = SM::sumPyrogenThresh;
    
    AgeGroup::init (monitoring);
    
    // by default, none are active
    active.reset ();
    SurveyMeasureMap codeMap;
    
    scnXml::OptionSet::OptionSequence sOSeq = monitoring.getSurveyOptions().getOption();
    for (scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
	active[codeMap[it->getName()]] = it->getValue();
    }
}
void AgeGroup::init (const scnXml::Monitoring& monitoring) {
    const scnXml::AgeGroup::GroupSequence& groups = monitoring.getAgeGroup().getGroup();
    /* note that the last age group includes individuals who are        *
    * either younger than Lowerbound or older than the last Upperbound */
    _upperbound.resize (groups.size() + 1);
    _lowerbound = monitoring.getAgeGroup().getLowerbound();
    if (!(_lowerbound <= 0.0))
	throw util::xml_scenario_error ("Expected survey age-group lowerbound of 0");

    for (size_t i = 0;i < groups.size(); ++i) {
	_upperbound[i] = groups[i].getUpperbound();
    }
    _upperbound[groups.size()] = numeric_limits<double>::infinity();
}

void AgeGroup::update (double ageYears) {
    while (ageYears > _upperbound[index]){
	++index;
    }
}


// -----  Non-static members  -----


Survey::Survey(){}

void Survey::allocate(){
    size_t numAgeGroups = AgeGroup::getNumGroups();
    m_humanReportsInt.resize(boost::extents[Report::MI_NUM][numAgeGroups]);
    m_humanReportsDouble.resize(boost::extents[Report::MD_NUM][numAgeGroups]);
    
    m_nTransmit = numeric_limits<double>::signaling_NaN();
    m_annAvgK = numeric_limits<double>::signaling_NaN();
    // These 3 are set as a whole and so don't require allocation:
//   _innoculationsPerDayOfYear;
//   _kappaPerDayOfYear;
//   _innoculationsPerAgeGroup;
    
    m_inputEIR =numeric_limits<double>::signaling_NaN() ;
    m_simulatedEIR = numeric_limits<double>::signaling_NaN();
    
    m_Clinical_RDTs = 0;
    m_Clinical_Microscopy = 0;
}

void Survey::writeSummaryArrays (ostream& outputFile, int survey)
{
    size_t nAgeGroups = m_humanReportsInt.shape()[1] - 1;   // Don't write out last age-group
    for( size_t intMeasure = 0; intMeasure < Report::MI_NUM; ++intMeasure ){
        if( active[intReportMappings[intMeasure]] ){
            for( size_t j = 0; j < nAgeGroups; ++j ){
                outputFile << survey << "\t" << j + 1 << "\t" << intReportMappings[intMeasure];
                outputFile << "\t" << m_humanReportsInt[intMeasure][j] << lineEnd;
            }
        }
    }
    for( size_t dblMeasure = 0; dblMeasure < Report::MD_NUM; ++dblMeasure ){
        if( active[dblReportMappings[dblMeasure]] ){
            for( size_t j = 0; j < nAgeGroups; ++j ){
                outputFile << survey << "\t" << j + 1 << "\t" << dblReportMappings[dblMeasure];
                outputFile << "\t" << m_humanReportsDouble[dblMeasure][j] << lineEnd;
            }
        }
    }
    
  if (active[SM::nTransmit]) {
    writeValue (outputFile, SM::nTransmit, survey, m_nTransmit);
  }
  if (active[SM::annAvgK]) {
    writeValue (outputFile, SM::annAvgK, survey, m_annAvgK);
  }

  if (active[SM::innoculationsPerAgeGroup]) {
    for (int j = 0; j < (int) m_inoculationsPerAgeGroup.size() - 1; ++j) { // Don't write out last age-group
        outputFile << survey << "\t" << j + 1 << "\t" << SM::innoculationsPerAgeGroup;
        outputFile << "\t" << m_inoculationsPerAgeGroup[j] << lineEnd;
    }
  }
  
  if (active[SM::Vector_Nv0]) {
    writeMap (outputFile, SM::Vector_Nv0, survey, data_Vector_Nv0);
  }
  if (active[SM::Vector_Nv]) {
    writeMap (outputFile, SM::Vector_Nv, survey, data_Vector_Nv);
  }
  if (active[SM::Vector_Ov]) {
    writeMap (outputFile, SM::Vector_Ov, survey, data_Vector_Ov);
  }
  if (active[SM::Vector_Sv]) {
    writeMap (outputFile, SM::Vector_Sv, survey, data_Vector_Sv);
  }
  if (active[SM::inputEIR]) {
    writeValue (outputFile, SM::inputEIR, survey, m_inputEIR);
  }
  if (active[SM::simulatedEIR]) {
    writeValue (outputFile, SM::simulatedEIR, survey, m_simulatedEIR);
  }
  if (active[SM::Clinical_RDTs]) {
      writeValue (outputFile, SM::Clinical_RDTs, survey, m_Clinical_RDTs);
  }
  if (active[SM::Clinical_DrugUsage]) {
      writeMap (outputFile, SM::Clinical_DrugUsage, survey, m_Clinical_DrugUsage);
  }
  if (active[SM::Clinical_DrugUsageIV]) {
      writeMap (outputFile, SM::Clinical_DrugUsageIV, survey, m_Clinical_DrugUsageIV);
  }
  if (active[SM::Clinical_Microscopy]) {
      writeValue (outputFile, SM::Clinical_Microscopy, survey, m_Clinical_Microscopy);
  }
}


template <class T>
void writeValue (ostream& file, int measure, int survey, T& value)
{
    file << survey << "\t" << 0 << "\t" << measure;
    file << "\t" << value << lineEnd;
}

void writeMap (ostream& file, int measure, int survey, map<string,double>& data)
{
    for (map<string,double>::const_iterator it = data.begin(); it != data.end(); ++it) {
        file << survey << "\t" << it->first << "\t" << measure;
        file << "\t" << it->second << lineEnd;
    }
}


void Survey::checkpoint(istream& stream){
    BOOST_STATIC_ASSERT( ReportsIntAgeT::dimensionality == 2 );
    size_t len0, len1;
    len0 & stream;
    len1 & stream;
    if( len0 != Report::MI_NUM || len1 != AgeGroup::getNumGroups() )
        throw util::checkpoint_error( "wrong survey data size" );
    m_humanReportsInt.resize(boost::extents[len0][len1]);
    for( size_t i = 0; i < len0; ++i ){
        for( size_t j = 0; j < len1; ++j ){
            m_humanReportsInt[i][j] & stream;
        }
    }
    BOOST_STATIC_ASSERT( ReportsDblAgeT::dimensionality == 2 );
    len0 & stream;
    len1 & stream;
    if( len0 != Report::MD_NUM || len1 != AgeGroup::getNumGroups() )
        throw util::checkpoint_error( "wrong survey data size" );
    m_humanReportsDouble.resize(boost::extents[len0][len1]);
    for( size_t i = 0; i < len0; ++i ){
        for( size_t j = 0; j < len1; ++j ){
            m_humanReportsDouble[i][j] & stream;
        }
    }
}

void Survey::checkpoint(ostream& stream) const{
    BOOST_STATIC_ASSERT( ReportsIntAgeT::dimensionality == 2 );
    size_t len0 = m_humanReportsInt.shape()[0], len1 = m_humanReportsInt.shape()[1];
    len0 & stream;
    len1 & stream;
    for( size_t i = 0; i < len0; ++i ){
        for( size_t j = 0; j < len1; ++j ){
            m_humanReportsInt[i][j] & stream;
        }
    }
    BOOST_STATIC_ASSERT( ReportsDblAgeT::dimensionality == 2 );
    len0 = m_humanReportsDouble.shape()[0], len1 = m_humanReportsDouble.shape()[1];
    len0 & stream;
    len1 & stream;
    for( size_t i = 0; i < len0; ++i ){
        for( size_t j = 0; j < len1; ++j ){
            m_humanReportsDouble[i][j] & stream;
        }
    }
}

} }