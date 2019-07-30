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
#include "mon/management.h"

#include <cstdlib>
#include <boost/xpressive/xpressive.hpp>
#include <iomanip>

using namespace boost::xpressive;

namespace OM {

// ———  SimTime stuff  ———

// Scenario constants
int SimData::interval;
size_t SimData::steps_per_year;
double SimData::years_per_step;

SimDate sim::s_start;
SimDate sim::s_end;

SimTime sim::s_max_human_age;

// Global variables
#ifndef NDEBUG
bool sim::in_update = false;
#endif
SimTime sim::s_t0;
SimTime sim::s_t1;

SimTime sim::s_interv;

using util::CommandLine;


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

/** Returns SimDate::never() when it doesn't recognise a date. Throws when it does
 * but encounters definite format errors. */
SimDate parseDate( const std::string& str ){
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
        return SimDate::origin()
            + SimTime::fromYearsI(year)
            + SimTime::roundToTSFromDays(monthStart[month-1] + day - 1);
    }else{
        return SimDate::never();
    }
}

SimDate readDate( const std::string& str, DefaultUnit defUnit ){
    SimDate date = parseDate( str );
    if( date != SimDate::never() ){
        if( date > sim::startDate() + SimTime::fromYearsI(500) ){
            cerr << "Warning: date is a long time after start date. Did you forget to set monitoring/startDate?" << endl;
        }else if( date < sim::startDate() ){
            throw util::format_error( "date of event is before the start of monitoring" );
        }
        return date;
    }else{
        if( CommandLine::option(CommandLine::DEPRECATION_WARNINGS) ){
            cerr << "Deprecation warning: time specified via duration \"" << str << "\" where a date could be used; recommended to use a date (e.g. 2011-12-20)" << endl;
        }
        return sim::startDate() + readDuration( str, defUnit );
    }
}

}       // end of UnitParse namespace


ostream& operator<<( ostream& stream, SimTime time ){
    if( time.inDays() % sim::DAYS_IN_YEAR == 0 ){
        stream << time.inYears() << 'y';
    } else {
        stream << time.inDays() << 'd';
    }
    return stream;
}

ostream& operator<<( ostream& stream, SimDate date ){
    int days = date.d;
    if (days < 0) {
        // Shouldn't happen; best still to print something
        stream << days << 'd';
    } else {
        int year = days / sim::DAYS_IN_YEAR;
        days -= year * sim::DAYS_IN_YEAR;
        
        int month = 0;
        while( days >= UnitParse::monthStart[month+1] ) ++month;
        days -= UnitParse::monthStart[month];
        assert( month < 12 && days < UnitParse::monthLen[month] );
        
        // Inconsistency in year vs month and day: see parseDate()
        stream << setfill('0') << setw(4) << year << '-'
                << setw(2) << (month+1) << '-' << (days+1)
                << setw(0) << setfill(' ');
    }
    return stream;
}

void sim::init( const scnXml::Scenario& scenario ){
    SimData::interval = scenario.getModel().getParameters().getInterval();
    SimData::steps_per_year = SimTime::oneYear().inSteps();
    SimData::years_per_step = 1.0 / SimData::steps_per_year;
    
    sim::s_max_human_age =
        SimTime::fromYearsD( scenario.getDemography().getMaximumAgeYrs() );
    
    sim::s_start = SimDate::origin();
    auto mon = scenario.getMonitoring();
    if( mon.getStartDate().present() ){
        try{
            // on failure, this throws or returns SimDate::never()
            sim::s_start = UnitParse::parseDate( mon.getStartDate().get() );
            if( sim::s_start == SimDate::never() ){
                throw util::format_error( "invalid format (expected YYYY-MM-DD)" );
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("monitoring/startDate: ").append(e.message()) );
        }
    }
    
    sim::s_interv = SimTime::never();    // large negative number
    
    sim::s_end = mon::readSurveyDates( mon );
}

}
