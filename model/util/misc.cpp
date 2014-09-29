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

namespace OM {

// ———  Unit parsing stuff  ———
namespace UnitParse {

int longToInt( long x ){
    int y = static_cast<int>( x );
    if( static_cast<long>( y ) != x ){
        throw util::format_error( "underflow/overflow" );
    }
    return y;
}

SimTime daysToTSExact( long days ){
    SimTime t = sim::fromDays( longToInt(days) );
    if( mod(t.inDays(), sim::oneTS().inDays()) != 0 ){
        throw util::format_error( "input in days is not a whole number of time steps" );
    }
    return t;
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
            return daysToTSExact( v );
        }else if( defUnit == STEPS ){
            return sim::fromTS( longToInt(v) );
        }else{
            throw SWITCH_DEFAULT_EXCEPTION;
        }
    }else if( len + 1 == str.size() ){
        // one extra character found; is this a unit?
        char u = str[len];
        if( u == 'd' || u == 'D' ){
            return daysToTSExact( v );
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
            return daysToTSExact( v );
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
            return daysToTSExact( v );
        }else if( u == 't' || u == 'T' ){
            return sim::fromTS( longToInt(v) );
        }
        // otherwise, fall through to below
    }
    throw util::format_error( string("bad format: '").append(str).append("' (try e.g. 1 or 2d or 3s or 4y)") );
}

SimTime readDate( const std::string& str, DefaultUnit defUnit ){
    throw util::unimplemented_exception("reading dates");
}

}
}
