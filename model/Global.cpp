/*

 This file is part of OpenMalaria.
 
 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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
#include "util/errors.hpp"

namespace OM {
    int Global::interval;
    int Global::intervalsPer5Days;
    size_t Global::intervalsPerYear;
    double Global::yearsPerInterval;
    int Global::maxAgeIntervals;
    int Global::lifespanInitIntervals;
    
    int Global::simulationTime;
    int Global::timeStep;
    
    void Global::init () {
	interval = InputData().getModel().getParameters().getInterval();
	if (Global::DAYS_IN_YEAR % interval !=  0) {
	    cerr << "Global::DAYS_IN_YEAR not a multiple of interval" << endl;
	    exit(-12);
	}
	intervalsPer5Days = 5/interval;
	intervalsPerYear = Global::DAYS_IN_YEAR/interval;
	yearsPerInterval = double(interval) / double(Global::DAYS_IN_YEAR);
	
	double maxAgeYears = InputData().getDemography().getMaximumAgeYrs();
	
	// Changed in r756 (2010-03-04): Did cast max-age-years to int before multiplying (minor effect on output).
	maxAgeIntervals = static_cast<int> (maxAgeYears * intervalsPerYear);
	// For schema 22: Make sure first day of intervention period is same
	// day-of-year of first day of simulation.
	lifespanInitIntervals = static_cast<int>(std::ceil(maxAgeYears)) * intervalsPerYear;
	// NOTE: This is also because both transmission models need at least
	// one year for initialization. This change is partially the opposite of
	// r756/r770, but not really the same.
	if( lifespanInitIntervals < (int)intervalsPerYear )
	    throw util::xml_scenario_error( "maximumAgeYrs must be positive" );
    }
}
