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

SimTime readShortDuration( const std::string& str, DefaultUnit defUnit ){
    const char * s = str.c_str();
    char * end = 0;
    
    long v = std::strtol( s, &end, 10 );
    size_t len = end - s;
    if( len == str.size() ){
        // no unit given; examine our policy:
        if( util::ModelOptions::option(util::REQUIRE_UNITS) || defUnit == NONE ){
            throw util::format_error( "unit required but not given (try e.g. 5d or 12t)" );
        }else if( defUnit == DAYS ){
            return sim::fromDays( v );
        }else if( defUnit == STEPS ){
            return sim::fromTS( v );
        }else{
            throw SWITCH_DEFAULT_EXCEPTION;
        }
    }else if( len + 1 == str.size() ){
        // one extra character found; is this a unit?
        char u = str[len];
        if( u == 'd' || u == 'D' ){
            return sim::fromDays( v );
        }else if( u == 't' || u == 'T' ){
            return sim::fromTS( v );
        }
        // otherwise, fall through to below
    }
    throw util::format_error( string("bad format: '").append(str).append("' (try e.g. 1 or 2d or 3s)") );
}

SimTime readDuration( const std::string& str, DefaultUnit defUnit ){
    throw util::unimplemented_exception("long durations (rounding behaviour?)");
}

SimTime readDate( const std::string& str, DefaultUnit defUnit ){
    throw util::unimplemented_exception("reading dates");
}

}
}
