/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#ifndef Hmod_ContinuousInterventions
#define Hmod_ContinuousInterventions

#include "Global.h"
#include <limits>
#include <map>

class InterventionsPerAge {
    public:
	InterventionsPerAge (int timeStep) :
	    ageTimeSteps(timeStep),
	    covVaccine(numeric_limits<double>::quiet_NaN()),
	    covITN(numeric_limits<double>::quiet_NaN()),
	    covIPTI(numeric_limits<double>::quiet_NaN())
	{}
	
	inline int getAgeTimeSteps () {
	    return ageTimeSteps;
	}
	
        void addVaccine(double cov) {
	    covVaccine = cov;
	}
	void addITN(double cov) {
	    covITN = cov;
	}
	void addIPTI(double cov) {
	    covIPTI = cov;
	}
	
    private:
	InterventionsPerAge () {}	// ctor only for use by map; an assignment takes place after creation
	friend InterventionsPerAge& std::map<int,InterventionsPerAge>::operator[](const int&);
	
	int ageTimeSteps;
	double covVaccine;
	double covITN;
	double covIPTI;
};

/** A continuous intervention scheduler.
 *
 * Groups all continuously deployed interventions together, deploying each when the
 * host reaches the target age. */
class ContinuousInterventions {
    public:
	ContinuousInterventions () :
	    next(0)
	{}
	
	void deploy (int ageTS);
	
	static void initParameters ();
	
    private:
	size_t next;	// index
	
	static vector<InterventionsPerAge> intervs;
};
#endif
