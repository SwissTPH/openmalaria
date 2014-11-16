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

#include "Monitoring/Surveys.h"
#include "Monitoring/AgeGroup.h"
#include "Host/Human.h"
#include "mon/management.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "schema/monitoring.h"
#include "schema/scenario.h"    // TODO: only for analysisNo

#include <boost/static_assert.hpp>
#include <stdexcept>

namespace OM { namespace Monitoring {
using WithinHost::diagnostics;
 
// -----  Utility forward declarations  -----

template <class T>
void writeValue(ostream& file, int measure, int survey, T& value);

void writeMap (ostream& file, int measure, int survey, map<string,double>& data);

size_t intReportMappings[Report::MI_NUM];

// -----  Static members  -----

vector<SimTime> AgeGroup::upperBound;
bitset<SM::NUM_SURVEY_OPTIONS> Survey::active;
const Diagnostic* Survey::m_diagnostic = 0;

class SurveyMeasureMap {
    // Lookup table to translate the strings used in the XML file to the internal enumerated values:
    map<string,SM::SurveyMeasure> codeMap;
    set<string> removedCodes;
    
public:
    SurveyMeasureMap () {
        codeMap["nHost"] = SM::BLANK;
        codeMap["nInfect"] = SM::BLANK;
        codeMap["nExpectd"] = SM::BLANK;
        codeMap["nPatent"] = SM::BLANK;
        codeMap["sumLogPyrogenThres"] = SM::BLANK;
        codeMap["sumlogDens"] = SM::BLANK;
        codeMap["totalInfs"] = SM::BLANK;
        codeMap["totalPatentInf"] = SM::BLANK;
        codeMap["sumPyrogenThresh"] = SM::BLANK;
        codeMap["nTreatments1"] = SM::BLANK;
        codeMap["nTreatments2"] = SM::BLANK;
        codeMap["nTreatments3"] = SM::BLANK;
        codeMap["nUncomp"] = SM::BLANK;
        codeMap["nSevere"] = SM::BLANK;
        codeMap["nSeq"] = SM::BLANK;
        codeMap["nHospitalDeaths"] = SM::BLANK;
        codeMap["nIndDeaths"] = SM::BLANK;
        codeMap["nDirDeaths"] = SM::BLANK;
        codeMap["nHospitalRecovs"] = SM::BLANK;
        codeMap["nHospitalSeqs"] = SM::BLANK;
        codeMap["nNMFever"] = SM::BLANK;
        codeMap["Clinical_FirstDayDeaths"] = SM::BLANK;
        codeMap["Clinical_HospitalFirstDayDeaths"] = SM::BLANK;
        codeMap["nNmfDeaths"] = SM::BLANK;
        codeMap["sumAge"] = SM::BLANK;
        
        codeMap["nTransmit"] = SM::nTransmit;
        codeMap["nEPIVaccinations"] = SM::nEPIVaccinations;
        codeMap["allCauseIMR"] = SM::allCauseIMR;
        codeMap["nMassVaccinations"] = SM::nMassVaccinations;
        codeMap["annAvgK"] = SM::annAvgK;
        codeMap["innoculationsPerAgeGroup"] = SM::innoculationsPerAgeGroup;
        codeMap["Vector_Nv0"] = SM::Vector_Nv0;
        codeMap["Vector_Nv"] = SM::Vector_Nv;
        codeMap["Vector_Ov"] = SM::Vector_Ov;
        codeMap["Vector_Sv"] = SM::Vector_Sv;
        codeMap["inputEIR"] = SM::inputEIR;
        codeMap["simulatedEIR"] = SM::simulatedEIR;
        codeMap["Clinical_RDTs"] = SM::Clinical_RDTs;
        codeMap["nNewInfections"] = SM::nNewInfections;
        codeMap["nMassITNs"] = SM::nMassITNs;
        codeMap["nEPI_ITNs"] = SM::nEPI_ITNs;
        codeMap["nMassIRS"] = SM::nMassIRS;
        codeMap["nMassGVI"] = SM::nMassGVI;
        codeMap["Clinical_Microscopy"] = SM::Clinical_Microscopy;
        codeMap["nMDAs"] = SM::nMDAs;
        codeMap["nMassScreenings"] = SM::nMassScreenings;
        codeMap["nCtsIRS"] = SM::nCtsIRS;
        codeMap["nCtsGVI"] = SM::nCtsGVI;
        codeMap["nCtsMDA"] = SM::nCtsMDA;
        codeMap["nCtsScreenings"] = SM::nCtsScreenings;
        codeMap["nSubPopRemovalTooOld"] = SM::nSubPopRemovalTooOld;
        codeMap["nSubPopRemovalFirstEvent"] = SM::nSubPopRemovalFirstEvent;
        codeMap["nPQTreatments"] = SM::nPQTreatments;
        codeMap["nTreatDiagnostics"] = SM::nTreatDiagnostics;
        codeMap["nMassRecruitOnly"] = SM::nMassRecruitOnly;
        codeMap["nCtsRecruitOnly"] = SM::nCtsRecruitOnly;
        codeMap["nTreatDeployments"] = SM::nTreatDeployments;
        
        removedCodes.insert("contrib");
        removedCodes.insert("nIPTDoses");
        removedCodes.insert("nAddedToCohort");
        removedCodes.insert("nRemovedFromCohort");
        removedCodes.insert("nAntibioticTreatments");
    }
    
    SM::SurveyMeasure operator[] (const string s) {
        map<string,SM::SurveyMeasure>::iterator codeIt = codeMap.find (s);
        if (codeIt == codeMap.end()) {
            ostringstream msg;
            if( removedCodes.count(s) > 0 ) msg << "Removed";
            else msg << "Unrecognised";
            msg << " survey option: \"" << s << '"';
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


void Survey::init( const OM::Parameters& parameters,
                   const scnXml::Scenario& scenario,
                   const scnXml::Monitoring& monitoring,
                   size_t nSurveys ){
    intReportMappings[Report::MI_VACCINATION_TIMED] = SM::nMassVaccinations;
    intReportMappings[Report::MI_VACCINATION_CTS] = SM::nEPIVaccinations;
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
    intReportMappings[Report::MI_N_SP_REM_TOO_OLD] = SM::nSubPopRemovalTooOld;
    intReportMappings[Report::MI_N_SP_REM_FIRST_EVENT] = SM::nSubPopRemovalFirstEvent;
    intReportMappings[Report::MI_PQ_TREATMENTS] = SM::nPQTreatments;
    intReportMappings[Report::MI_TREAT_DIAGNOSTICS] = SM::nTreatDiagnostics;
    intReportMappings[Report::MI_RECRUIT_TIMED] = SM::nMassRecruitOnly;
    intReportMappings[Report::MI_RECRUIT_CTS] = SM::nCtsRecruitOnly;
    intReportMappings[Report::MI_TREAT_DEPLOYMENTS] = SM::nTreatDeployments;
    
    AgeGroup::init (monitoring);
    
    mon::initialise( nSurveys, AgeGroup::getNumGroups(), Surveys.numCohortSets(), monitoring );
    
    // by default, none are active
    active.reset ();
    SurveyMeasureMap codeMap;
    
    scnXml::OptionSet::OptionSequence sOSeq = monitoring.getSurveyOptions().getOption();
    for (scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
	active[codeMap[it->getName()]] = it->getValue();
    }
    
    const scnXml::Surveys& surveys = monitoring.getSurveys();
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        // So far the implemented Vivax code does not produce parasite
        // densities, thus this diagnostic model cannot be used.
        m_diagnostic = &diagnostics::make_deterministic( numeric_limits<double>::quiet_NaN() );
    }else if( surveys.getDetectionLimit().present() ){
        if( surveys.getDiagnostic().present() ){
            throw util::xml_scenario_error( "monitoring/surveys: do not "
                "specify both detectionLimit and diagnostic" );
        }
        if( util::CommandLine::option( util::CommandLine::DEPRECATION_WARNINGS ) ){
            std::cerr << "Deprecation warning: monitoring/surveys: "
                "specification of \"diagnostic\" is suggested over \"detectionLimit\"" << std::endl;
        }
        
        // This controls whether the detection limit is specified relative to
        // the Garki or other methods.
        double densitybias = numeric_limits<double>::quiet_NaN();
        if (util::ModelOptions::option (util::GARKI_DENSITY_BIAS)) {
            densitybias = parameters[Parameters::DENSITY_BIAS_GARKI];
        } else {
            if( scenario.getAnalysisNo().present() ){
                int analysisNo = scenario.getAnalysisNo().get();
                if ((analysisNo >= 22) && (analysisNo <= 30)) {
                    cerr << "Warning: these analysis numbers used to mean "
                        "use Garki density bias. If you do want to use this, "
                        "specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong." << endl;
                }
            }
            densitybias = parameters[Parameters::DENSITY_BIAS_NON_GARKI];
        }
        double detectionLimit = surveys.getDetectionLimit().get() * densitybias;
        m_diagnostic = &diagnostics::make_deterministic( detectionLimit );
    }else{
        if( !surveys.getDiagnostic().present() ){
            throw util::xml_scenario_error( "monitoring/surveys: require "
                "either detectionLimit or diagnostic" );
        }
        if( util::ModelOptions::option(util::GARKI_DENSITY_BIAS) ){
            throw util::xml_scenario_error( "Use of GARKI_DENSITY_BIAS is not "
                "appropriate when monitoring/surveys/diagnostic is used." );
        }
        m_diagnostic = &diagnostics::get( surveys.getDiagnostic().get() );
    }
}
void AgeGroup::init (const scnXml::Monitoring& monitoring) {
    const scnXml::MonAgeGroup::GroupSequence& groups = monitoring.getAgeGroup().getGroup();
    if (!(monitoring.getAgeGroup().getLowerbound() <= 0.0))
        throw util::xml_scenario_error ("Expected survey age-group lowerbound of 0");
    
    // The last age group includes individuals too old for reporting
    upperBound.resize( groups.size() + 1 );
    for (size_t i = 0;i < groups.size(); ++i) {
        // convert to SimTime, rounding down to the next time step
        upperBound[i] = sim::fromYearsD( groups[i].getUpperbound() );
    }
    upperBound[groups.size()] = sim::future();
}

void AgeGroup::update (SimTime age) {
    while (age >= upperBound[index]){
        ++index;
    }
}


// -----  Non-static members  -----


Survey::Survey() :
        m_nTransmit(numeric_limits<double>::signaling_NaN()),
        m_annAvgK(numeric_limits<double>::signaling_NaN()),
        m_inputEIR(numeric_limits<double>::signaling_NaN()),
        m_simulatedEIR(numeric_limits<double>::signaling_NaN()),
        m_Clinical_RDTs(numeric_limits<int>::min()),
        m_Clinical_Microscopy(numeric_limits<int>::min())
{}

void Survey::allocate(){
    size_t numAgeGroups = AgeGroup::getNumGroups();
    size_t nCohortSets = Surveys.numCohortSets();
    m_humanReportsInt.resize(boost::extents[Report::MI_NUM][numAgeGroups][nCohortSets]);
    
    
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

Survey& Survey::addInt( ReportMeasureI measure, const Host::Human &human, int val ){
    size_t ageIndex = human.getMonitoringAgeGroup().i();
#ifndef NDEBUG
    if( static_cast<size_t>(measure.code) >= m_humanReportsInt.shape()[0] ||
        ageIndex >= m_humanReportsInt.shape()[1] ||
        human.cohortSet() >= m_humanReportsInt.shape()[2]
    ){
        cout << "Index out of bounds for survey:\t" << static_cast<void*>(this)
            << "\nmeasure number\t" << measure.code << " of " << m_humanReportsInt.shape()[0]
            << "\nage group\t" << ageIndex << " of " << m_humanReportsInt.shape()[1]
            << "\ncohort set\t" << human.cohortSet() << " of " << m_humanReportsInt.shape()[2]
            << endl;
    }
#endif
    m_humanReportsInt[measure.code][ageIndex][human.cohortSet()] += val;
    return *this;
}

Survey& Survey::addInt( ReportMeasureI measure, AgeGroup ageGroup, uint32_t cohortSet, int val ){
    if( static_cast<size_t>(measure.code) >= m_humanReportsInt.shape()[0] ||
        ageGroup.i() >= m_humanReportsInt.shape()[1] ){
        cout << "Index out of bounds:\n"
            "survey\t" << static_cast<void*>(this)
            << "\nalloc\t" << m_humanReportsInt.shape()[0] << "\t" << m_humanReportsInt.shape()[1]
            << "\nindex\t" << measure.code << "\t" << ageGroup.i() << endl;
    }
    m_humanReportsInt[measure.code][ageGroup.i()][cohortSet] += val;
    return *this;
}

void Survey::writeSummaryArrays (ostream& outputFile, int survey)
{
    size_t nAgeGroups = m_humanReportsInt.shape()[1] - 1;   // Don't write out last age-group
    size_t nCohortSets = m_humanReportsInt.shape()[2];
    mon::writeMHI( outputFile, survey );
    for( size_t intMeasure = 0; intMeasure < Report::MI_NUM; ++intMeasure ){
        if( active[intReportMappings[intMeasure]] ){
            for( size_t cohortSet = 0; cohortSet < nCohortSets; ++cohortSet ){
                for( size_t ageGroup = 0; ageGroup < nAgeGroups; ++ageGroup ){
                    // Yeah, >999 age groups clashes with cohort sets, but unlikely a real issue
                    size_t ageCohortId = 1000 * Surveys.cohortSetOutputId( cohortSet ) + ageGroup + 1;
                    outputFile << survey << "\t" << ageCohortId << "\t" << intReportMappings[intMeasure];
                    outputFile << "\t" << m_humanReportsInt[intMeasure][ageGroup][cohortSet] << lineEnd;
                }
            }
        }
    }
    mon::writeMHD( outputFile, survey );
    
  if (active[SM::nTransmit]) {
    writeValue (outputFile, SM::nTransmit, survey, m_nTransmit);
  }
  if (active[SM::annAvgK]) {
    writeValue (outputFile, SM::annAvgK, survey, m_annAvgK);
  }

  if (active[SM::innoculationsPerAgeGroup]) {
    for (int ageGroup = 0; ageGroup < (int) m_inoculationsPerAgeGroup.size() - 1; ++ageGroup) { // Don't write out last age-group
        outputFile << survey << "\t" << ageGroup + 1 << "\t" << SM::innoculationsPerAgeGroup;
        outputFile << "\t" << m_inoculationsPerAgeGroup[ageGroup] << lineEnd;
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
    BOOST_STATIC_ASSERT( ReportsIntAgeT::dimensionality == 3 );
    size_t len0, len1, len2;
    len0 & stream;
    len1 & stream;
    len2 & stream;
    if( len0 != Report::MI_NUM || len1 != AgeGroup::getNumGroups() ||
        len2 != Surveys.numCohortSets() )
    {
        throw util::checkpoint_error( "wrong survey data size" );
    }
    m_humanReportsInt.resize(boost::extents[len0][len1][len2]);
    for( size_t i = 0; i < len0; ++i ){
        for( size_t j = 0; j < len1; ++j ){
            for( size_t k = 0; k < len2; ++k ){
                m_humanReportsInt[i][j][k] & stream;
            }
        }
    }
}

void Survey::checkpoint(ostream& stream) const{
    BOOST_STATIC_ASSERT( ReportsIntAgeT::dimensionality == 3 );
    size_t len0 = m_humanReportsInt.shape()[0],
        len1 = m_humanReportsInt.shape()[1],
        len2 = m_humanReportsInt.shape()[2];
    len0 & stream;
    len1 & stream;
    len2 & stream;
    for( size_t i = 0; i < len0; ++i ){
        for( size_t j = 0; j < len1; ++j ){
            for( size_t k = 0; k < len2; ++k ){
                m_humanReportsInt[i][j][k] & stream;
            }
        }
    }
}

} }