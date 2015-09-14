/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

#include "mon/info.h"
#include "mon/management.h"
#include "mon/AgeGroup.h"
#include "mon/reporting.h"
#include "interventions/InterventionManager.hpp"
#include "util/BoincWrapper.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "util/timeConversions.h"
#include "schema/monitoring.h"

#include "WithinHost/Diagnostic.h"

#include <gzstream/gzstream.h>
#include <fstream>
#include <boost/algorithm/string.hpp>

namespace OM {
namespace mon {

// ———  surveys  ———

namespace impl{
    // Constants or defined during init:
    size_t nSurveys = 0;
    size_t nCohortSets = 1;     // default: just the whole population
    vector<SimTime> surveyTimes;        // times of all surveys (from XML)
}

void initSurveyTimes( const OM::Parameters& parameters,
                   const scnXml::Scenario& scenario,
                   const scnXml::Monitoring& monitoring ){
    const scnXml::Surveys::SurveyTimeSequence& survs =
        monitoring.getSurveys().getSurveyTime();
    
    impl::surveyTimes.reserve (survs.size());   // insufficient reservation if repetitions are used
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
                    impl::surveyTimes.push_back( cur );
                    cur += step;
                }
            }else{
                impl::surveyTimes.push_back( cur );
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("surveyTime: ").append(e.message()) );
        }
    }
    // sort times:
    sort( impl::surveyTimes.begin(), impl::surveyTimes.end() );
    // remove duplicates:
    vector<SimTime>::iterator newEnd = unique( impl::surveyTimes.begin(), impl::surveyTimes.end() );
    impl::nSurveys = distance(impl::surveyTimes.begin(), newEnd);
    if( impl::nSurveys < impl::surveyTimes.size() ){
        std::cerr << "Warning: " << (impl::surveyTimes.size() - impl::nSurveys)
                << " duplicate survey times omitted. Survey numbers do not "
                "include these." << std::endl
                << "Note: the OpenMalaria v33 release will not work correctly "
                "with this XML." << std::endl;
    }
    impl::surveyTimes.resize( impl::nSurveys );
    
    if( util::CommandLine::option( util::CommandLine::PRINT_SURVEY_TIMES ) ){
        bool haveDate = UnitParse::haveDate();
        std::cout << "Survey\tsteps\tdays";
        if( haveDate ) std::cout << "\tdate";
        std::cout << std::endl;
        for( size_t i = 0; i < impl::nSurveys; ++i ){
            std::cout << (i+1) << '\t' << impl::surveyTimes[i].inSteps()
                    << '\t' << impl::surveyTimes[i].inDays();
            if( haveDate ) std::cout << '\t' << impl::surveyTimes[i];
            std::cout << std::endl;
        }
    }
    
    if( monitoring.getCohorts().present() ){
        // this needs to be set early, but we can't set cohortSubPopIds until after InterventionManager is initialised
        impl::nCohortSets = static_cast<uint32_t>(1) << monitoring.getCohorts().get().getSubPop().size();
    }
    
    mon::AgeGroup::init( monitoring );

    internal::initReporting( scenario );
}

void initMainSim(){
    if( impl::nSurveys > 0 ){
        impl::currentSurvey = 0;
        impl::nextSurveyTime = impl::surveyTimes[impl::currentSurvey];
    }
}
void concludeSurvey(){
    impl::currentSurvey += 1;
    if( impl::currentSurvey >= impl::nSurveys ){
        // After the last survey has completed:
        impl::currentSurvey = NOT_USED;
        impl::nextSurveyTime = sim::future();
    }else{
        impl::nextSurveyTime = impl::surveyTimes[impl::currentSurvey];
    }
}

SimTime nextSurveyTime(){
    return impl::nextSurveyTime;
}
SimTime finalSurveyTime(){
    return impl::surveyTimes[impl::surveyTimes.size()-1];
}

void writeSurveyData ()
{
#ifdef WITHOUT_BOINC
    ofstream outputFile;          // without boinc, use plain text (for easy reading)
#else
    ogzstream outputFile;         // with, use gzip
#endif
    
    // This locale ensures uniform formatting of nans and infs on all platforms.
    std::locale old_locale;
    std::locale nfn_put_locale(old_locale, new boost::math::nonfinite_num_put<char>);
    outputFile.imbue( nfn_put_locale );

    string output_filename = util::BoincWrapper::resolveFile(
        util::CommandLine::getOutputName() );
    
    outputFile.open( output_filename.c_str(), std::ios::out | std::ios::binary );
    
    outputFile.width (0);
    // For additional control:
    // outputFile.precision (6);
    // outputFile << scientific;
    
    internal::write( outputFile );
    
    outputFile.close();
}


// ———  AgeGroup  ———

vector<SimTime> AgeGroup::upperBound;

void AgeGroup::init (const scnXml::Monitoring& monitoring) {
    const scnXml::MonAgeGroup::GroupSequence& groups =
        monitoring.getAgeGroup().getGroup();
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


// ———  Cohort sets  ———

using interventions::ComponentId;
vector<uint32_t> cohortSubPopNumbers;   // value is output number
map<ComponentId,uint32_t> cohortSubPopIds;      // value is internal index (used above)

bool notPowerOfTwo( uint32_t num ){
    for( uint32_t i = 0; i <= 21; ++i ){
        if( num == (static_cast<uint32_t>(1) << i) )
            return false;
    }
    return true;
}
// Init cohort sets. Depends on interventions (initialise those first).
void initCohorts( const scnXml::Monitoring& monitoring )
{
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

uint32_t updateCohortSet( uint32_t old, ComponentId subPop, bool isMember ){
    map<ComponentId,uint32_t>::const_iterator it = cohortSubPopIds.find( subPop );
    if( it == cohortSubPopIds.end() ) return old;       // sub-pop not used in cohorts
    uint32_t subPopId = static_cast<uint32_t>(1) << it->second;        // 1 bit positive
    return (old & ~subPopId) | (isMember ? subPopId : 0);
}

uint32_t internal::cohortSetOutputId(uint32_t cohortSet){
    uint32_t outNum = 0;
    assert( (cohortSet >> cohortSubPopNumbers.size()) == 0 );
    for( uint32_t i = 0; i < cohortSubPopNumbers.size(); ++i ){
        if( cohortSet & (static_cast<uint32_t>(1) << i) ){
            outNum += cohortSubPopNumbers[i];
        }
    }
    return outNum;
}

} }
