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

#include "Global.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/UnitParse.h"
#include "schema/scenario.h"
#include "mon/management.h"

#include <cstdlib>
#include <iomanip>
#include <regex>

namespace OM {

// ———  SimTime stuff  ———

// Scenario constants
int SimData::interval;
size_t SimData::steps_per_year;
double SimData::years_per_step;

SimTime sim::s_start;
SimTime sim::s_end;

SimTime sim::s_max_human_age;

// Global variables
#ifndef NDEBUG
bool sim::in_update = false;
#endif
SimTime sim::s_t0;
SimTime sim::s_t1;

SimTime sim::s_interv;

using util::CommandLine;

ostream& operator<<( ostream& stream, SimTime time ){
    if( time.inDays() % sim::DAYS_IN_YEAR == 0 ){
        stream << time.inYears() << 'y';
    } else {
        stream << time.inDays() << 'd';
    }
    return stream;
}

// ostream& operator<<( ostream& stream, SimTime date ){
//     int days = date.d;
//     if (days < 0) {
//         // Shouldn't happen; best still to print something
//         stream << days << 'd';
//     } else {
//         int year = days / sim::DAYS_IN_YEAR;
//         days -= year * sim::DAYS_IN_YEAR;
        
//         int month = 0;
//         while( days >= UnitParse::monthStart[month+1] ) ++month;
//         days -= UnitParse::monthStart[month];
//         assert( month < 12 && days < UnitParse::monthLen[month] );
        
//         // Inconsistency in year vs month and day: see parseDate()
//         stream << setfill('0') << setw(4) << year << '-'
//                 << setw(2) << (month+1) << '-' << (days+1)
//                 << setw(0) << setfill(' ');
//     }
//     return stream;
// }

void sim::init( const scnXml::Scenario& scenario ){
    SimData::interval = scenario.getModel().getParameters().getInterval();
    SimData::steps_per_year = sim::oneYear().inSteps();
    SimData::years_per_step = 1.0 / SimData::steps_per_year;
    
    sim::s_max_human_age =
        sim::fromYearsD( scenario.getDemography().getMaximumAgeYrs() );
    
    sim::s_start = sim::origin();
    auto mon = scenario.getMonitoring();
    if( mon.getStartDate().present() ){
        try{
            // on failure, this throws or returns sim::never()
            sim::s_start = UnitParse::parseDate( mon.getStartDate().get() );
            if( sim::s_start == sim::never() ){
                throw util::format_error( "invalid format (expected YYYY-MM-DD)" );
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("monitoring/startDate: ").append(e.message()) );
        }
    }
    
    sim::s_interv = sim::never();    // large negative number
    
    sim::s_end = mon::readSurveyDates( mon );
}

}
