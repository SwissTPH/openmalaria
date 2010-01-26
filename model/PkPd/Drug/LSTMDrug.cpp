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
	double IC50 = typeData->parameters.IC50;
	double slope = typeData->parameters.slope;
	double elimination_rate_constant = typeData->parameters.elimination_rate_constant;
	double max_killing_rate = typeData->parameters.max_killing_rate;
	double vol_dist = typeData->parameters.vol_dist;
	
	
	double totalFactor = 1.0;   /* KW-	The drug factor being passed to melissa - this begins with a value of 1, it assumes no drug affect is seen
										this vaule is updated in the for loop, value decreases with increasing drug effect. */
    double startTime = 0.0;		/* KW-	Use the information from medicate to determine the time elapsed from 0 to first dose.
										Use the information on dose timings from medicate to update this value at the end of the "for" loop . */
    for (list<Dose>::iterator dose = doses.begin(); dose!=doses.end(); ++dose) {
	double duration = dose->time - startTime;
	 
	//KW - start concentration is equal to the end concentration of the previous time step
	double conc_after_decay = concentration * exp(elimination_rate_constant *  duration);			//KW - need to find time.
	double drug_effect = pow( (pow(IC50,slope) + pow(conc_after_decay,slope)) / (pow(IC50,slope) + pow(concentration,slope)), max_killing_rate / (elimination_rate_constant * slope));
	totalFactor *= drug_effect;	
	startTime = dose->time;		// KW - Increment the time 

	concentration = conc_after_decay + dose->mg / (vol_dist*weight_kg);
								// KW - The concentration value here is moved to the top of this loop and used in the conc_feter_decay and drug_effect equations
    }
   
	double duration = 24*60 - startTime;
	 
	//KW - start concentration is equal to the end concentration of the previous time step
	double conc_after_decay = concentration * exp(elimination_rate_constant *  duration);			//KW - need to find time.
	double drug_effect = pow( (pow(IC50,slope) + pow(conc_after_decay,slope)) / (pow(IC50,slope) + pow(concentration,slope)), max_killing_rate / (elimination_rate_constant * slope));
	totalFactor *= drug_effect;	

	concentration = conc_after_decay;
	
	doses.clear ();				// KW - Clear doses to ensure they don't interfer with those given on the next day.
    
    return totalFactor;			/* KW -	Returning drug effect per day, per drug
										Where is this returned to? 
										Need to multiply the return values of all drugs in one day together to get one totalFactor value */
}

double LSTMDrug::decayFactor (double time) {
    //TBD
    return 0.0;	// TODO (best always return _something_, even if nonsense)
}

} }