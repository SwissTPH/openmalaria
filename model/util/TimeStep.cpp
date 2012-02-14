/*
 This file is part of OpenMalaria.

 Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <cstdlib>
#include <cmath>

#include "Global.h"
#include "inputData.h"
#include "util/errors.h"

namespace OM {
namespace util {

// Days in a year. Should be a multiple of interval.
const int DAYS_IN_YEAR = 365;

TimeStep::TimeStep() : _ts(-0x3FFFFFFF) {}
TimeStep TimeStep::fromDaysNearest( double d ){
    return TimeStep( static_cast<int>(std::floor( d / interval + 0.5 )) );
}
TimeStep TimeStep::fromDays( double d ){
    return TimeStep( static_cast<int>(std::floor( d / interval )) );
}

TimeStep::ReadOnly<int> TimeStep::interval;
TimeStep::ReadOnly<double> TimeStep::yearsPerInterval;
TimeStep TimeStep::intervalsPer5Days( TimeStep::never );
TimeStep TimeStep::intervalsPerYear( TimeStep::never );
TimeStep TimeStep::maxAgeIntervals( TimeStep::never );
TimeStep::ReadOnly<int> TimeStep::stepsPerYear;

const TimeStep TimeStep::never, TimeStep::future(0x7FFFFFFF);

TimeStep TimeStep::simulation( TimeStep::never );
TimeStep TimeStep::interventionPeriod( TimeStep::never );

void TimeStep::init( int daysPerTimeStep, double maxAgeYears ) {
    interval = daysPerTimeStep;
    if (DAYS_IN_YEAR % interval != 0) {
        throw util::xml_scenario_error( "Global::DAYS_IN_YEAR not a multiple of interval" );
    }else if (maxAgeYears < 1.0){
        throw util::xml_scenario_error( "maximumAgeYrs must be at least 1" );
    }
    intervalsPer5Days = TimeStep(5/interval);
    stepsPerYear = DAYS_IN_YEAR/interval;
    intervalsPerYear = TimeStep( stepsPerYear );
    yearsPerInterval = double(interval) / DAYS_IN_YEAR;
    
    // Maximum age of an individual:
    maxAgeIntervals = TimeStep( static_cast<int>(maxAgeYears * stepsPerYear) );
}

}
}

