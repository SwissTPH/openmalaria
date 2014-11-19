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

#include "mon/info.h"
#include "mon/management.h"
#include "mon/AgeGroup.h"
#include "mon/reporting.h"
#include "Monitoring/Survey.h"  // for init
#include "util/errors.h"
#include "schema/monitoring.h"
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
    
    impl::surveyTimes.reserve (survs.size());
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
    sort( impl::surveyTimes.begin(), impl::surveyTimes.end() );
    impl::nSurveys = impl::surveyTimes.size();
    
    if( monitoring.getCohorts().present() ){
        // this needs to be set early, but we can't set cohortSubPopIds until after InterventionManager is initialised
        impl::nCohortSets = static_cast<uint32_t>(1) << monitoring.getCohorts().get().getSubPop().size();
    }
    
    Monitoring::Survey::init( parameters, scenario, monitoring );
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

} }
