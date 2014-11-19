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

// Must not be included _after_ boost/math/special_functions/fpclassify.hpp
#include <boost/math/nonfinite_num_facets.hpp>

#include "Monitoring/Surveys.h"
#include "Simulator.h"
#include "Clinical/CaseManagementCommon.h"
#include "interventions/InterventionManager.hpp"
#include "mon/management.h"
#include "util/BoincWrapper.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "schema/monitoring.h"

#include <gzstream/gzstream.h>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace OM { namespace Monitoring {
using interventions::ComponentId;

SurveysType Surveys;
size_t Survey::m_surveyNumber;
Survey* Survey::m_current;
uint32_t nCohortSets = 1;     // default: just the whole population
vector<uint32_t> cohortSubPopNumbers;   // value is output number
map<ComponentId,uint32_t> cohortSubPopIds;      // value is internal index (used above)

bool notPowerOfTwo( uint32_t num ){
    for( uint32_t i = 0; i <= 21; ++i ){
        if( num == (static_cast<uint32_t>(1) << i) )
            return false;
    }
    return true;
}
void SurveysType::init( const OM::Parameters& parameters,
                   const scnXml::Scenario& scenario,
                   const scnXml::Monitoring& monitoring ){
    Survey::m_surveyNumber = 0;
    if( monitoring.getCohorts().present() ){
        // this needs to be set early, but we can't set cohortSubPopIds until after InterventionManager is initialised
        nCohortSets = static_cast<uint32_t>(1) << monitoring.getCohorts().get().getSubPop().size();
    }
    
    const scnXml::Surveys::SurveyTimeSequence& survs = monitoring.getSurveys().getSurveyTime();

    m_surveysTimeIntervals.reserve (survs.size() + 1);
    for (size_t i = 0; i < survs.size(); ++i) {
        const scnXml::SurveyTime& surv = survs[i];
        try{
            std::string s = surv;
            boost::algorithm::trim(s);
            SimTime cur = UnitParse::readDate( s, UnitParse::STEPS );
            if( surv.getRepeatStep().present() != surv.getRepeatEnd().present() ){
                throw util::xml_scenario_error( "surveyTime: use of repeatStep or repeatEnd without other" );
            }
            if( surv.getRepeatStep().present() ){
                SimTime step = UnitParse::readDuration( surv.getRepeatStep().get(), UnitParse::NONE );
                if( step < sim::oneTS() ){
                    throw util::xml_scenario_error( "surveyTime: repeatStep must be >= 1" );
                }
                SimTime end = UnitParse::readDate( surv.getRepeatEnd().get(), UnitParse::NONE );
                while(cur < end){
                    m_surveysTimeIntervals.push_back( cur );
                    cur += step;
                }
            }else{
                m_surveysTimeIntervals.push_back( cur );
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("surveyTime: ").append(e.message()) );
        }
    }
    sort( m_surveysTimeIntervals.begin(), m_surveysTimeIntervals.end() );
    m_surveysTimeIntervals.push_back( sim::never() );
    m_nextSurveyTime = m_surveysTimeIntervals[0];

    Survey::init( parameters, scenario, monitoring );

    m_surveys.resize (m_surveysTimeIntervals.size());
    if( !Simulator::isCheckpoint() ){
        for (size_t i = 0; i < m_surveys.size(); ++i)
            m_surveys[i].allocate();
    }
    Survey::m_current = &m_surveys[0];
}

void SurveysType::init2( const scnXml::Monitoring& monitoring, size_t nSpecies )
{
    mon::initialise( m_surveysTimeIntervals.size()-1, Surveys.numCohortSets(),
                     nSpecies, monitoring );
    
    if( monitoring.getCohorts().present() ){
        const scnXml::Cohorts monCohorts = monitoring.getCohorts().get();
        uint32_t nextId = 0;
        for( scnXml::Cohorts::SubPopConstIterator it = monCohorts.getSubPop().begin(),
            end = monCohorts.getSubPop().end(); it != end; ++it )
        {
            ComponentId compId = interventions::InterventionManager::getComponentId( it->getId() );
            bool inserted = cohortSubPopIds.insert( make_pair(compId,nextId) ).second;
            if( !inserted ){
                throw util::xml_scenario_error(
                    string("cohort specification uses sub-population \"").append(it->getId())
                    .append("\" more than once") );
            }
            if( it->getNumber() < 0 || notPowerOfTwo( it->getNumber() ) ){
                throw util::xml_scenario_error(
                    string( "cohort specification assigns sub-population \"").append(it->getId())
                    .append("\" a number which is not a power of 2 (up to 2^21)") );
            }
            cohortSubPopNumbers.push_back( it->getNumber() );
            nextId += 1;
        }
    }
}

Survey& Survey::getSurvey(size_t n){
    assert (n < Surveys.m_surveys.size());
    return Surveys.m_surveys[n];
}
SimTime Survey::getLastSurveyTime () {
    return Surveys.m_surveysTimeIntervals[Surveys.m_surveys.size()-2];   // last entry in this list is sim::never()
}

uint32_t SurveysType::numCohortSets()const{
    return nCohortSets;
}

uint32_t Survey::updateCohortSet( uint32_t old, ComponentId subPop, bool isMember ){
    map<ComponentId,uint32_t>::const_iterator it = cohortSubPopIds.find( subPop );
    if( it == cohortSubPopIds.end() ) return old;       // sub-pop not used in cohorts
    uint32_t subPopId = static_cast<uint32_t>(1) << it->second;        // 1 bit positive
    return (old & ~subPopId) | (isMember ? subPopId : 0);
}

uint32_t SurveysType::cohortSetOutputId(uint32_t cohortSet) const{
    uint32_t outNum = 0;
    assert( (cohortSet >> cohortSubPopNumbers.size()) == 0 );
    for( uint32_t i = 0; i < cohortSubPopNumbers.size(); ++i ){
        if( cohortSet & (static_cast<uint32_t>(1) << i) ){
            outNum += cohortSubPopNumbers[i];
        }
    }
    return outNum;
}

void SurveysType::incrementSurveyPeriod()
{
  m_nextSurveyTime = m_surveysTimeIntervals[Survey::m_surveyNumber];
  ++Survey::m_surveyNumber;
  if (Survey::m_surveyNumber >= m_surveys.size())
    // In this case, m_nextSurveyTime gets set to sim::never() so no further surveys get taken
    Survey::m_surveyNumber = 0;
  Survey::m_current = &m_surveys[Survey::m_surveyNumber];
}

void SurveysType::writeSummaryArrays ()
{
#ifdef WITHOUT_BOINC
  ofstream outputFile;		// without boinc, use plain text (for easy reading)
#else
  ogzstream outputFile;		// with, use gzip
#endif
    
    // This locale ensures uniform formatting of nans and infs on all platforms.
    locale old_locale;
    locale nfn_put_locale(old_locale, new boost::math::nonfinite_num_put<char>);
    outputFile.imbue( nfn_put_locale );

    string output_filename = util::BoincWrapper::resolveFile (util::CommandLine::getOutputName());

  outputFile.open (output_filename.c_str(), ios::out | ios::binary);

  outputFile.width (0);
  // For additional control:
  //   outputFile.precision (6);
  //   outputFile << scientific;

  for (size_t i = 1; i < m_surveys.size(); ++i){
    mon::write1( outputFile, i );
    m_surveys[i].writeSummaryArrays (outputFile, i);
    mon::write2( outputFile, i );
  }

  //Infant mortality rate is a single number, therefore treated separately
  // Note: Storing a single value instead of one per reporting period is inconsistent with other
  // reporting, but I believe required for parameterisation.
  if (Survey::active[SM::allCauseIMR]) {
    outputFile << 1 << "\t" << 1 << "\t" << SM::allCauseIMR;
    outputFile << "\t" << Clinical::infantAllCauseMort() << lineEnd;
  }

  outputFile.close();
}

void SurveysType::checkpoint (istream& stream) {
    m_nextSurveyTime & stream;
    m_surveysTimeIntervals & stream;
    Survey::m_surveyNumber & stream;
    // read those surveys checkpointed, call allocate on the rest:
    for( size_t i = 1; i <= Survey::m_surveyNumber; ++i )
        m_surveys[i] & stream;
    m_surveys[0].allocate();
    for( size_t i = Survey::m_surveyNumber + 1; i < m_surveys.size(); ++i )
        m_surveys[i].allocate();
    Survey::m_current = &m_surveys[Survey::m_surveyNumber];
}
void SurveysType::checkpoint (ostream& stream) {
    m_nextSurveyTime & stream;
    m_surveysTimeIntervals & stream;
    Survey::m_surveyNumber & stream;
    // checkpoint only those surveys used; exclude 0 since that's a "write only DB"
    for( size_t i = 1; i <= Survey::m_surveyNumber; ++i )
        m_surveys[i] & stream;
}

} }
