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
#include "util/timeConversions.h"
#include "schema/scenario.h"

#include <cstdlib>
#include <boost/xpressive/xpressive.hpp>
#include <iomanip>

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

using util::CommandLine;


TimeDisplayHelper SimTime::date() const {
    return TimeDisplayHelper(*this, TimeDisplayHelper::DATE);
}
/// Display as a duration with automatic units
TimeDisplayHelper SimTime::duration() const {
    return TimeDisplayHelper(*this, TimeDisplayHelper::DURATION);
}

// ———  Unit parsing stuff  ———
namespace UnitParse {

static int monthLen[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
static int monthStart[] = { // offsets of start of month from start of year
    0 /*Jan*/, 31, 59, 90,
    120 /*May*/, 151, 181, 212,
    243 /*Sept*/, 273, 304, 334,
    365 /*next year: stop condition in formatDate*/
};

template<typename T>
int castToInt( T x ){
    int y = static_cast<int>( x );
    if( static_cast<T>( y ) != x ){
        throw util::format_error( "underflow/overflow" );
    }
    return y;
}

bool haveDate(){
    return interv_start_date != SimTime::never();
}

SimTime readShortDuration( const std::string& str, DefaultUnit defUnit ){
    const char * s = str.c_str();
    char * end = 0;
    
    long v = std::strtol( s, &end, 10 );
    size_t len = end - s;
    if( len == str.size() ){
        // no unit given; examine our policy:
        if( v == 0 ){   // don't need a unit in this case
            return SimTime::zero();
        }else if( defUnit == NONE ){
            throw util::format_error( "unit required but not given (try e.g. 5d or 12t)" );
        }else{
            if( CommandLine::option(CommandLine::DEPRECATION_WARNINGS) ){
                cerr << "Deprecation warning: duration \"" << str << "\" specified without unit; it is recommended to do so (e.g. 5d or 1t)" << endl;
            }
            if( defUnit == DAYS ){
                return SimTime::roundToTSFromDays( v );
            }else if( defUnit == STEPS ){
                return SimTime::fromTS( castToInt(v) );
            }else{
                throw SWITCH_DEFAULT_EXCEPTION;
            }
        }
    }else if( len + 1 == str.size() ){
        // one extra character found; is this a unit?
        char u = str[len];
        if( u == 'd' || u == 'D' ){
            return SimTime::roundToTSFromDays( v );
        }else if( u == 't' || u == 'T' ){
            return SimTime::fromTS( castToInt(v) );
        }
        // otherwise, fall through to below
    }
    throw util::format_error( string("bad format: '").append(str).append("' (try e.g. 1 or 2d or 3s)") );
}

// Parses str, returns v and may set unit
double parseDurationAndUnit( const std::string& str, DefaultUnit& unit ){
    const char * s = str.c_str();
    char * end = 0;
    
    double v = std::strtod( s, &end );
    size_t len = end - s;
    if( len + 1 == str.size() ){
        // one extra character found; is this a unit?
        char u = str[len];
        if( u == 'y' || u == 'Y' ) unit = YEARS;
        else if( u == 'd' || u == 'D' ) unit = DAYS;
        else if( u == 't' || u == 'T' ) unit = STEPS;
        else {
            throw util::format_error( string("unknown unit: '").append(str).append("' (try e.g. 1 or 2d or 3s or 4y)") );
        }
    }else if( len == str.size() ){
        if( v == 0.0 ) unit = DAYS;        // special case: 0 does not require a unit; pretend the unit is days
        else if( unit == NONE /*no default set for this value*/ ){
            throw util::format_error( "unit required but not given (try e.g. 5d or 12t or 2.3y)" );
        }else if( CommandLine::option(CommandLine::DEPRECATION_WARNINGS) ){
            cerr << "Deprecation warning: duration \"" << str << "\" specified without unit; it is recommended to do so (e.g. 5d or 1t or 0.5y)" << endl;
        }
    }else{
        throw util::format_error( string("bad format: '").append(str).append("' (try e.g. 1 or 2d or 3s or 4y)") );
    }
    return v;
}

SimTime readDuration( const std::string& str, DefaultUnit unit ){
    double v = parseDurationAndUnit( str, unit );
    
    if( unit == YEARS ) return SimTime::fromYearsN( v );
    else if( v != std::floor(v) ){
        throw util::format_error( "fractional values are only allowed when the unit is years (e.g. 0.25y)" );
    }else if( unit == DAYS ) return SimTime::roundToTSFromDays( v );
	else if( unit == STEPS ) return SimTime::fromTS( castToInt(v) );
    else throw SWITCH_DEFAULT_EXCEPTION;
}

double durationToDays( const std::string& str, DefaultUnit unit ){
    double v = parseDurationAndUnit( str, unit );
    
    if( unit == YEARS ) return v * SimTime::oneYear().inDays();
    else if( unit == DAYS ) return v;
    else if( unit == STEPS ) return v * SimTime::oneTS().inDays();
    else throw SWITCH_DEFAULT_EXCEPTION;
}

/** Returns SimTime::never() when it doesn't recognise a date. Throws when it does
 * but encounters definite format errors. */
SimTime parseDate( const std::string& str ){
    static sregex dateRx = sregex::compile ( "(\\d{4})-(\\d{1,2})-(\\d{1,2})" );
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
        // Inconsistency: time "zero" is 0000-01-01, not 0001-01-01. Since
        // dates are always relative to another date, the extra year doesn't
        // actually affect anything.
        return SimTime::fromYearsI(year) + SimTime::roundToTSFromDays(monthStart[month-1] + day - 1);
    }else{
        return SimTime::never();
    }
}

SimTime readDate( const std::string& str, DefaultUnit defUnit ){
    SimTime date = parseDate( str );
    if( date != SimTime::never() ){
        if( interv_start_date == SimTime::never() ){
            throw util::format_error( "must set monitoring/startDate when using dates" );
        }else if( date < interv_start_date ){
            throw util::format_error( "date of event is before the start of monitoring" );
        }else {
            return date - interv_start_date;  // externally, intervention dates are relative to the start
        }
    }else{
        if( CommandLine::option(CommandLine::DEPRECATION_WARNINGS) ){
            cerr << "Deprecation warning: time specified via duration \"" << str << "\" where a date could be used; recommended to use a date (e.g. 2011-12-20)" << endl;
        }
        return readDuration( str, defUnit );
    }
}

}       // end of UnitParse namespace


ostream& operator<<( ostream& stream, TimeDisplayHelper timeDisplay ){
    if( timeDisplay.mode == TimeDisplayHelper::DATE && interv_start_date != SimTime::never() ){
        // format as a date
        SimTime date = interv_start_date + timeDisplay.time;
        assert( date >= SimTime::zero() );
        int year = date / SimTime::oneYear();
        assert( year >= 0 );
        int remainder = (date - SimTime::oneYear() * year).inDays();
        int month = 0;
        while( remainder >= UnitParse::monthStart[month+1] ) ++month;
        remainder -= UnitParse::monthStart[month];
        assert( month < 12 && remainder < UnitParse::monthLen[month] );
        // Inconsistency in year vs month and day: see parseDate()
        stream << setfill('0') << setw(4) << year << '-'
                << setw(2) << (month+1) << '-' << (remainder+1)
                << setw(0) << setfill(' ');
    } else {
        // format as a duration
        if( timeDisplay.time.inDays() % SimTime::DAYS_IN_YEAR == 0 ){
            stream << timeDisplay.time.inYears() << 'y';
        } else {
            stream << timeDisplay.time.inDays() << 'd';
        }
    }
    return stream;
}

void sim::init( const scnXml::Scenario& scenario ){
    SimTime::interval = scenario.getModel().getParameters().getInterval();
    SimTime::steps_per_year = SimTime::oneYear().inSteps();
    SimTime::years_per_step = 1.0 / SimTime::steps_per_year;
    sim::max_human_age = SimTime::fromYearsD( scenario.getDemography().getMaximumAgeYrs() );
    if( scenario.getMonitoring().getStartDate().present() ){
        try{
            // on failure, this throws or returns SimTime::never()
            interv_start_date = UnitParse::parseDate( scenario.getMonitoring().getStartDate().get() );
            if( interv_start_date == SimTime::never() ){
                throw util::format_error( "invalid format (expected YYYY-MM-DD)" );
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("monitoring/startDate: ").append(e.message()) );
        }
    }
    
    sim::interv_time = SimTime::never();    // large negative number
}

}
