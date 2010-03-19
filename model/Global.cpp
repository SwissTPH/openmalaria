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

namespace OM {
    int Global::interval;
    int Global::intervalsPer5Days;
    size_t Global::intervalsPerYear;
    double Global::yearsPerInterval;
    int Global::maxAgeIntervals;
    
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
	// Changed in r756 (2010-03-04): Did cast max-age-years to int before multiplying (minor effect on output).
	maxAgeIntervals = static_cast<int> (InputData().getMaximumAgeYrs() * intervalsPerYear);
    }
}
