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

#include "PkPd/Drug/LSTMDrug.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrug::LSTMDrug(const LSTMDrugType* type) : Drug(),
    typeData (type),
    concentration (0.0)
{}


void LSTMDrug::storeDose (double concentration, int delay) {
    doses.push_back (Dose (concentration, delay));
}


double LSTMDrug::calculateDrugFactor(const ProteomeInstance* infProteome, double ageYears, double weight_kg) {
    typeData->parameters.max_killing_rate;
    
    double totalFactor = 1.0;
    double startTime = 0.0;
    for (list<Dose>::iterator i = doses.begin(); i!=doses.end(); ++i) {
	//TODO::
	/* calculate t (from time startTime to this drug time)
	 calculate kill factor over this t
	 */
	totalFactor *= 1.0 /*FACTOR*/;
	
	startTime = i->time;
	concentration += i->mg;
    }
    //TODO: and here (from startTime to end time)
    doses.clear ();
    
    return totalFactor;
}

double LSTMDrug::decayFactor (double time) {
    //TBD
    return 0.0;	// TODO (best always return _something_, even if nonsense)
}

} }