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
    
LSTMDrug::LSTMDrug(const LSTMDrugType* type) :
    typeData (type),
    concentration (0.0)
{}


void LSTMDrug::storeDose (double time, double qty) {
    doses.push_back (Dose (time, qty));
}


/* Function to avoid repeating some operations in calculateDrugFactor().
 * @param duration Time-span acted over with units days.
 * @returns a survival factor (no units). */
inline double drugEffect (const LSTMDrugPDParameters& PD_params, double elimination_rate_constant, double& concentration, double duration) {
    //KW - start concentration is equal to the end concentration of the previous time step
    double conc_after_decay = concentration * exp(elimination_rate_constant *  duration);
    
    // Note: these look a little different from original equations because PD_params.IC50_pow_slope
    // and PD_params.power are calculated when read from the scenario document instead of here.
    double numerator = PD_params.IC50_pow_slope + pow(conc_after_decay,PD_params.slope);
    double denominator = PD_params.IC50_pow_slope + pow(concentration,PD_params.slope);
    double drug_effect = pow( numerator / denominator, PD_params.power );
    
    concentration = conc_after_decay;
    return drug_effect;
}

double LSTMDrug::calculateDrugFactor(uint32_t proteome_ID, double ageYears, double weight_kg) const {
    //FIXME: check units of time
    
    double totalFactor = 1.0;		/* KW-	The drug factor being passed to melissa - this begins with a value of 1, it assumes no drug affect is seen
									this vaule is updated in the for loop, value decreases with increasing drug effect. */
    double startTime = 0.0;		/* KW-	Use the information from medicate to determine the time elapsed from 0 to first dose.
									Use the information on dose timings from medicate to update this value at the end of the for loop.
									Run drugEffect function after for loop to find drug effect from last dose to end of day. */
    
    double dist_weight_inv = 1.0 / (typeData->vol_dist * weight_kg);
    // Make a copy of concetration and use that over today. Don't adjust concentration because this
    // function may be called multiple times (or not at all) in a day.
    double concentration_today = concentration;
    
    size_t allele = (proteome_ID >> typeData->allele_rshift) & typeData->allele_mask;
    const LSTMDrugPDParameters& PD_params = typeData->PD_params[allele];
    
    for (list<Dose>::const_iterator dose = doses.begin(); dose!=doses.end(); ++dose) {
	double duration = dose->time - startTime;
	totalFactor *= drugEffect (PD_params, typeData->elimination_rate_constant, concentration_today, duration);	
	concentration_today += dose->mg * dist_weight_inv;
	
	startTime = dose->time;		// KW - Increment the time (assuming doses are in order of time)
    }
   
    double duration = 1.0 - startTime;
    totalFactor *= drugEffect (PD_params, typeData->elimination_rate_constant, concentration_today, duration);	
    
    return totalFactor;			/* KW -	Returning drug effect per day, per drug */
}

bool LSTMDrug::updateConcentration (double weight_kg) {
    //FIXME: check units of time
    
    double startTime = 0.0;		// as in calculateDrugFactor()
    double dist_weight_inv = 1.0 / (typeData->vol_dist * weight_kg);
    
    for (list<Dose>::const_iterator dose = doses.begin(); dose!=doses.end(); ++dose) {
	double duration = dose->time - startTime;
	concentration *= exp(typeData->elimination_rate_constant *  duration);
	concentration += dose->mg * dist_weight_inv;
	
	startTime = dose->time;	
    }
   
    double duration = 1.0 - startTime;
    concentration *= exp(typeData->elimination_rate_constant *  duration);
    
    doses.clear ();				// Clear today's dose list — they've been added to concentration now.
    // return true when concentration is no longer significant:
    return concentration < typeData->negligible_concentration;
}

} }