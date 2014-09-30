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

#include "Global.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "schema/scenario.h"

#include <cstdlib>
#include <boost/xpressive/xpressive.hpp>

using namespace boost::xpressive;

namespace OM {

// ———  SimTime stuff  ———

// global constants:
int SimTime::interval;
size_t SimTime::steps_per_year;
double SimTime::years_per_step;
SimTime sim::max_human_age;
SimTime interv_start_date;

// global variables:
#ifndef NDEBUG
bool sim::in_update = false;
#endif
SimTime sim::time0;
SimTime sim::time1;
SimTime sim::interv_time;


// ———  Unit parsing stuff  ———
namespace UnitParse {

int longToInt( long x ){
    int y = static_cast<int>( x );
    if( static_cast<long>( y ) != x ){
        throw util::format_error( "underflow/overflow" );
    }
    return y;
}

SimTime readShortDuration( const std::string& str, DefaultUnit defUnit ){
    const char * s = str.c_str();
    char * end = 0;
    
    long v = std::strtol( s, &end, 10 );
    size_t len = end - s;
    if( len == str.size() ){
        // no unit given; examine our policy:
        if( v == 0 ){   // don't need a unit in this case
            return sim::zero();
        }else if( util::ModelOptions::option(util::REQUIRE_UNITS) || defUnit == NONE ){
            throw util::format_error( "unit required but not given (try e.g. 5d or 12t)" );
        }else if( defUnit == DAYS ){
            return sim::roundToTSFromDays( v );
        }else if( defUnit == STEPS ){
            return sim::fromTS( longToInt(v) );
        }else{
            throw SWITCH_DEFAULT_EXCEPTION;
        }
    }else if( len + 1 == str.size() ){
        // one extra character found; is this a unit?
        char u = str[len];
        if( u == 'd' || u == 'D' ){
            return sim::roundToTSFromDays( v );
        }else if( u == 't' || u == 'T' ){
            return sim::fromTS( longToInt(v) );
        }
        // otherwise, fall through to below
    }
    throw util::format_error( string("bad format: '").append(str).append("' (try e.g. 1 or 2d or 3s)") );
}

SimTime readDuration( const std::string& str, DefaultUnit defUnit ){
    const char * s = str.c_str();
    char * end = 0;
    
    double v = std::strtod( s, &end );
    size_t len = end - s;
    if( len == str.size() ){
        // no unit given; examine our policy:
        if( v == 0 ){   // don't need a unit in this case
            return sim::zero();
        }else if( util::ModelOptions::option(util::REQUIRE_UNITS) || defUnit == NONE ){
            throw util::format_error( "unit required but not given (try e.g. 5d or 12t or 2.3y)" );
        }else if( defUnit == YEARS ){
            return sim::fromYearsN( v );
        }else if( v != std::floor(v) ){
            throw util::format_error( "fractional values are only allowed when the unit is years (e.g. 0.25y)" );
        }else if( defUnit == DAYS ){
            return sim::roundToTSFromDays( v );
        }else if( defUnit == STEPS ){
            return sim::fromTS( longToInt(v) );
        }else{
            throw SWITCH_DEFAULT_EXCEPTION;
        }
    }else if( len + 1 == str.size() ){
        // one extra character found; is this a unit?
        char u = str[len];
        if( u == 'y' || u == 'Y' ){
            return sim::fromYearsN( v );
        }else if( v != std::floor(v) ){
            throw util::format_error( "fractional values are only allowed when the unit is years (e.g. 0.25y)" );
        }else if( u == 'd' || u == 'D' ){
            return sim::roundToTSFromDays( v );
        }else if( u == 't' || u == 'T' ){
            return sim::fromTS( longToInt(v) );
        }
        // otherwise, fall through to below
    }
    throw util::format_error( string("bad format: '").append(str).append("' (try e.g. 1 or 2d or 3s or 4y)") );
}

/** Returns sim::never() when it doesn't recognise a date. Throws when it does
 * but encounters definite format errors. */
SimTime parseDate( const std::string& str ){
    static sregex dateRx = sregex::compile ( "(\\d{4})-(\\d{1,2})-(\\d{1,2})" );
    static int monthLen[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    static int monthStart[] = { // offsets of start of month from start of year
        0 /*Jan*/, 31, 59, 90,
        120 /*May*/, 151, 181, 212,
        243 /*Sept*/, 273, 304, 334
    };
    smatch rx_what;
    if( regex_match(str, rx_what, dateRx) ){
        int year = std::atoi(rx_what[1].str().c_str()),
            month = std::atoi(rx_what[2].str().c_str()),
            day = std::atoi(rx_what[3].str().c_str());
        if( month < 1 || month > 12 || day < 1 || day > monthLen[month-1] ){
            throw util::format_error( string(str).append(" does not look like a date (expected YYYY-MM-DD with 1≤MM≤12 and 1≤DD≤(days in month))") );
        }
        // We can't overflow: the largest year possible is 9999 which is about
        // 3.6 million days; an int is good to more than 1e9.
        // We should also round to the nearest step
        return sim::fromYearsI(year) + sim::roundToTSFromDays(monthStart[month-1] + day - 1);
    }else{
        return sim::never();
    }
}

SimTime readDate( const std::string& str, DefaultUnit defUnit ){
    SimTime date = parseDate( str );
    if( date != sim::never() ){
        if( interv_start_date == sim::never() ){
            throw util::format_error( "must set monitoring/startDate when using dates" );
        }else if( date < interv_start_date ){
            throw util::format_error( "date of event is before the start of monitoring" );
        }else {
            return date - interv_start_date;  // externally, intervention dates are relative to the start
        }
    }else{
        if( util::ModelOptions::option(util::REQUIRE_DATES) ){
            throw util::format_error( string("expected to find a date, not ").append(str) );
        }else{
            return readDuration( str, defUnit );
        }
    }
}

}       // end of UnitParse namespace


void sim::init( const scnXml::Scenario& scenario ){
    SimTime::interval = scenario.getModel().getParameters().getInterval();
    SimTime::steps_per_year = sim::oneYear().inSteps();
    SimTime::years_per_step = 1.0 / SimTime::steps_per_year;
    sim::max_human_age = sim::fromTS( scenario.getDemography().getMaximumAgeYrs() * SimTime::steps_per_year );
    if( scenario.getMonitoring().getStartDate().present() ){
        try{
            // on failure, this throws or returns sim::never()
            interv_start_date = UnitParse::parseDate( scenario.getMonitoring().getStartDate().get() );
            if( interv_start_date == sim::never() ){
                throw util::format_error( "invalid format (expected YYYY-MM-DD)" );
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("monitoring/startDate: ").append(e.message()) );
        }
    }
    
    sim::interv_time = sim::never();    // large negative number
}

}
