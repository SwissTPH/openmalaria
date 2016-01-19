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

struct SurveyTime {
    SimTime time;       // time of survey
    size_t num; // if NOT_USED, the survey is not reported; if greater, this is the survey number
    
    /// Construct
    SurveyTime(SimTime time, size_t num) : time(time), num(num) {}
    
    inline bool isReported() const { return num != NOT_USED; }
};

namespace impl{
    // Constants or defined during init:
    size_t nSurveys = 0;        // number of reported surveys
    size_t nCohortSets = 1;     // default: just the whole population
    extern size_t surveyIndex;     // index in surveyTimes of next survey
    vector<SurveyTime> surveyTimes;     // times of surveys
}

void updateConditions();        // defined in mon.cpp

void initSurveyTimes( const OM::Parameters& parameters,
                   const scnXml::Scenario& scenario,
                   const scnXml::Monitoring& monitoring ){
    const scnXml::Surveys::SurveyTimeSequence& survs =
        monitoring.getSurveys().getSurveyTime();
    
    map<SimTime, bool> surveys;        // times of all surveys (from XML) and whether these are reporting
    
    for(size_t i = 0; i < survs.size(); ++i) {
        const scnXml::SurveyTime& surv = survs[i];
        try{
            std::string s = surv;
            boost::algorithm::trim(s);
            SimTime cur = UnitParse::readDate( s, UnitParse::STEPS );
            bool reporting = surv.getReported();
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
                    if( reporting ){
                        surveys[cur] = true;
                    }else{
                        surveys.insert(make_pair(cur, false));        // does not override existing pair with key 'cur'
                    }
                    cur += step;
                }
            }else{
                if( reporting ){
                    surveys[cur] = true;
                }else{
                    surveys.insert(make_pair(cur, false));        // does not override existing pair with key 'cur'
                }
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("surveyTime: ").append(e.message()) );
        }
    }
    
    impl::surveyTimes.clear();
    impl::surveyTimes.reserve(surveys.size());
    size_t n = 0;
    for( map<SimTime, bool>::const_iterator it = surveys.begin(); it != surveys.end(); ++it ){
        size_t num = NOT_USED;
        if( it->second ){
            num = n;
            n += 1;
        }
        impl::surveyTimes.push_back(SurveyTime(it->first, num));
    }
    impl::nSurveys = n;
    
    if( impl::surveyTimes.size() > 0 && !impl::surveyTimes.back().isReported() ){
        std::cerr << "Warning: the last survey is unreported. Having surveys beyond the last reported survey is pointless." << std::endl;
    }
    
    if( util::CommandLine::option( util::CommandLine::PRINT_SURVEY_TIMES ) ){
        bool haveDate = UnitParse::haveDate();
        std::cout << "Survey\tsteps\tdays";
        if( haveDate ) std::cout << "\tdate";
        std::cout << std::endl;
        for( size_t i = 0; i < impl::surveyTimes.size(); ++i ){
            const SurveyTime& survTime = impl::surveyTimes[i];
            if( !survTime.isReported() ) continue;
            std::cout << (survTime.num+1) << '\t' << survTime.time.inSteps()
                    << '\t' << survTime.time.inDays();
            if( haveDate ) std::cout << '\t' << survTime.time;
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

void updateSurveyNumbers() {
    if( impl::surveyIndex >= impl::surveyTimes.size() ){
        impl::survNumEvent = NOT_USED;
        impl::survNumStat = NOT_USED;
        impl::nextSurveyTime = sim::future();
    }else{
        for( size_t i = impl::surveyIndex; i < impl::surveyTimes.size(); ++i ){
            impl::survNumEvent = impl::surveyTimes[i].num;  // set to survey number or NOT_USED; this happens at least once!
            if( impl::survNumEvent != NOT_USED ) break;        // stop at first reported survey
        }
        const SurveyTime& nextSurvey = impl::surveyTimes[impl::surveyIndex];
        impl::survNumStat = nextSurvey.num;     // may be NOT_USED; this is intended
        impl::nextSurveyTime = nextSurvey.time;
    }
}
void initMainSim(){
    impl::surveyIndex = 0;
    impl::isInit = true;
    updateSurveyNumbers();
}
void concludeSurvey(){
    updateConditions();
    impl::surveyIndex += 1;
    updateSurveyNumbers();
}

SimTime nextSurveyTime(){
    return impl::nextSurveyTime;
}
SimTime finalSurveyTime(){
    return impl::surveyTimes[impl::surveyTimes.size()-1].time;
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
    for(size_t i = 0;i < groups.size(); ++i) {
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
