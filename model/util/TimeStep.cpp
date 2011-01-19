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
    return TimeStep( std::floor( d / interval + 0.5 ) );
}
TimeStep TimeStep::fromDays( double d ){
    return TimeStep( std::floor( d / interval ) );
}

TimeStep::ReadOnly<int> TimeStep::interval;
TimeStep::ReadOnly<double> TimeStep::yearsPerInterval;
TimeStep TimeStep::intervalsPer5Days( TimeStep::never );
TimeStep TimeStep::intervalsPerYear( TimeStep::never );
TimeStep TimeStep::maxAgeIntervals( TimeStep::never );
TimeStep TimeStep::lifespanInitIntervals( TimeStep::never );
TimeStep::ReadOnly<int> TimeStep::stepsPerYear;

const TimeStep TimeStep::never, TimeStep::future(0x7FFFFFFF);

TimeStep TimeStep::simulation( TimeStep::never );
TimeStep TimeStep::interventionPeriod( TimeStep::never );

void TimeStep::init( int daysPerTimeStep, double maxAgeYears ) {
    interval = daysPerTimeStep;
    if (DAYS_IN_YEAR % interval != 0) {
        cerr << "Global::DAYS_IN_YEAR not a multiple of interval" << endl;
        exit(-12);
    }
    intervalsPer5Days = TimeStep(5/interval);
    stepsPerYear = DAYS_IN_YEAR/interval;
    intervalsPerYear = TimeStep( stepsPerYear );
    yearsPerInterval = double(interval) / DAYS_IN_YEAR;

    // Changed in r756 (2010-03-04): Did cast max-age-years to int before multiplying (minor effect on output).
    maxAgeIntervals = TimeStep (maxAgeYears * stepsPerYear);
    // For schema 22: Make sure first day of intervention period is same
    // day-of-year of first day of simulation.
    lifespanInitIntervals = TimeStep(std::ceil(maxAgeYears) * stepsPerYear);
    // NOTE: This is also because both transmission models need at least
    // one year for initialization. This change is partially the opposite of
    // r756/r770, but not really the same.
    if ( lifespanInitIntervals < intervalsPerYear )
        throw util::xml_scenario_error( "maximumAgeYrs must be positive" );
}

}
}

