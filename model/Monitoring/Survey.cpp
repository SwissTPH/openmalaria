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
        codeMap["nEPIVaccinations"] = SM::BLANK;
        codeMap["nMassVaccinations"] = SM::BLANK;
        codeMap["nMassITNs"] = SM::BLANK;
        codeMap["nEPI_ITNs"] = SM::BLANK;
        codeMap["nMassIRS"] = SM::BLANK;
        codeMap["nMassGVI"] = SM::BLANK;
        codeMap["nMDAs"] = SM::BLANK;
        codeMap["nMassScreenings"] = SM::BLANK;
        codeMap["nCtsIRS"] = SM::BLANK;
        codeMap["nCtsGVI"] = SM::BLANK;
        codeMap["nCtsMDA"] = SM::BLANK;
        codeMap["nCtsScreenings"] = SM::BLANK;
        codeMap["nMassRecruitOnly"] = SM::BLANK;
        codeMap["nCtsRecruitOnly"] = SM::BLANK;
        codeMap["nTreatDeployments"] = SM::BLANK;
        codeMap["nNewInfections"] = SM::BLANK;
        codeMap["nSubPopRemovalTooOld"] = SM::BLANK;
        codeMap["nSubPopRemovalFirstEvent"] = SM::BLANK;
        codeMap["nPQTreatments"] = SM::BLANK;
        codeMap["nTreatDiagnostics"] = SM::BLANK;
        codeMap["nTransmit"] = SM::BLANK;
        codeMap["annAvgK"] = SM::BLANK;
        codeMap["inputEIR"] = SM::BLANK;
        codeMap["simulatedEIR"] = SM::BLANK;
        
        codeMap["allCauseIMR"] = SM::allCauseIMR;
        codeMap["innoculationsPerAgeGroup"] = SM::innoculationsPerAgeGroup;
        codeMap["Vector_Nv0"] = SM::Vector_Nv0;
        codeMap["Vector_Nv"] = SM::Vector_Nv;
        codeMap["Vector_Ov"] = SM::Vector_Ov;
        codeMap["Vector_Sv"] = SM::Vector_Sv;
        
        removedCodes.insert("Clinical_RDTs");
        removedCodes.insert("Clinical_Microscopy");
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


Survey::Survey() {}

void Survey::allocate(){
    // These 3 are set as a whole and so don't require allocation:
//   _innoculationsPerDayOfYear;
//   _kappaPerDayOfYear;
//   _innoculationsPerAgeGroup;
}

void Survey::writeSummaryArrays (ostream& outputFile, int survey)
{
    mon::write( outputFile, survey );
    
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

} }