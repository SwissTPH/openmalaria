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


inline double drugEffect (const LSTMDrugType& drugType, double& concentration, double duration, double weight_kg, double dose_mg) {
    //KW - start concentration is equal to the end concentration of the previous time step
    double conc_after_decay = concentration * exp(drugType.elimination_rate_constant *  duration);
    size_t allele = 0;	//FIXME: fudge using first allele until we get allele info.
    const LSTMDrugPDParameters& PD_params = drugType.PD_params[allele];
    
    double numerator = pow(PD_params.IC50,PD_params.slope) + pow(conc_after_decay,PD_params.slope);
    double denominator = pow(PD_params.IC50,PD_params.slope) + pow(concentration,PD_params.slope);
    double power = PD_params.max_killing_rate / (drugType.elimination_rate_constant * PD_params.slope);
    double drug_effect = pow( numerator / denominator, power );
    
    concentration = conc_after_decay + dose_mg / (drugType.vol_dist*weight_kg);
    
    return drug_effect;
}

double LSTMDrug::calculateDrugFactor(const ProteomeInstance* infProteome, double ageYears, double weight_kg) {
    double totalFactor = 1.0;		/* KW-	The drug factor being passed to melissa - this begins with a value of 1, it assumes no drug affect is seen
									this vaule is updated in the for loop, value decreases with increasing drug effect. */
    double startTime = 0.0;		/* KW-	Use the information from medicate to determine the time elapsed from 0 to first dose.
									Use the information on dose timings from medicate to update this value at the end of the for loop.
									Run drugEffect function after for loop to find drug effect from last dose to end of day. */
    for (list<Dose>::iterator dose = doses.begin(); dose!=doses.end(); ++dose) {
	double duration = dose->time - startTime;
	
	totalFactor *= drugEffect (*typeData, concentration, duration, weight_kg, dose->mg);	
	startTime = dose->time;		// KW - Increment the time 
    }
   
    double duration = 24*60 - startTime;
    
    totalFactor *= drugEffect (*typeData, concentration, duration, weight_kg, 0.0);	
    
    doses.clear ();				// KW - Clear doses to ensure they don't interfer with those given on the next day.
    
    return totalFactor;			/* KW -	Returning drug effect per day, per drug */
}

double LSTMDrug::decayFactor (double time) {
    //TBD
    return 0.0;	// TODO (best always return _something_, even if nonsense)
}

} }