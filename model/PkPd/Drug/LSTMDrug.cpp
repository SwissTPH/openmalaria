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
#include "util/errors.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrug::LSTMDrug(const LSTMDrugType& type) :
    typeData (&type),
    concentration (0.0)
{}


void LSTMDrug::medicate (double time, double qty, double weight) {
    double conc = qty / (typeData->vol_dist * weight);
    // multimap insertion: is ordered
    doses.insert (doses.end(), make_pair (time, conc));
}

void LSTMDrug::medicateIV (double duration, double endTime, double qty, double weight) {
    if( typeData->IV_params == NULL ){
	throw util::xml_scenario_error( "IV medication of a drug without IV parameters!" );
    }
    
    //TODO: decide whether we input mg/kg or just mg
    double mgPerKg = qty / weight;
    
    double infusRate = mgPerKg / duration;	// mg/day
    double elim_rate_const = typeData->IV_params->elimination_rate_constant;
    
    IV_doses.push_back( IV_dose( infusRate, duration ) );
    
    // TODO: check formula is correct:
    double conc = infusRate / (typeData->IV_params->vol_dist * elim_rate_const);
    conc *= (1.0 - exp( -elim_rate_const * duration ));
    
    // multimap insertion: is ordered
    doses.insert (doses.end(), make_pair (endTime, conc));
}


/* Function to avoid repeating some operations in calculateDrugFactor().
 * @param duration Time-span acted over with units days.
 * @returns a survival factor (no units). */
inline double drugEffect (const LSTMDrugPDParameters& PD_params, double neg_elimination_rate_constant, double& concentration, double duration) {
    //KW - start concentration is equal to the end concentration of the previous time step
    double conc_after_decay = concentration * exp(neg_elimination_rate_constant *  duration);
    
    // Note: these look a little different from original equations because PD_params.IC50_pow_slope
    // and PD_params.power are calculated when read from the scenario document instead of here.
    double numerator = PD_params.IC50_pow_slope + pow(conc_after_decay,PD_params.slope);
    double denominator = PD_params.IC50_pow_slope + pow(concentration,PD_params.slope);
    double drug_effect = pow( numerator / denominator, PD_params.power );
    
    concentration = conc_after_decay;
    return drug_effect;
}

// TODO: in high transmission, is this going to get called more than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrug::calculateDrugFactor(uint32_t proteome_ID) const {
    double totalFactor = 1.0;		/* KW-	The drug factor being passed to melissa - this begins with a value of 1, it assumes no drug effect is seen
									this vaule is updated in the for loop, value decreases with increasing drug effect. */
    double startTime = 0.0;			/* KW-	Use the information from medicate to determine the time elapsed from 0 to first dose.
									Use the information on dose timings from medicate to update this value at the end of the for loop.
									Run drugEffect function after for loop to find drug effect from last dose to end of day. */
    
    // Make a copy of concetration and use that over today. Don't adjust concentration because this
    // function may be called multiple times (or not at all) in a day.
    double concentration_today = concentration;
    
    size_t allele = (proteome_ID >> typeData->allele_rshift) & typeData->allele_mask;
    const LSTMDrugPDParameters& PD_params = typeData->PD_params[allele];
    
    for (multimap<double,double>::const_iterator dose = doses.begin(); dose!=doses.end(); ++dose) {
	if( dose->first >= 1.0 )
	    break;	// we know this and any more doses happen tomorrow; don't calculate factors now
	
	double duration = dose->first - startTime;
	// TODO: skip if duration == 0.0?
	totalFactor *= drugEffect (PD_params, typeData->neg_elimination_rate_constant, concentration_today, duration);
	concentration_today += dose->second;
	
	startTime = dose->first;		// KW - Increment the time (assuming doses are in order of time)
    }
   
    double duration = 1.0 - startTime;
    totalFactor *= drugEffect (PD_params, typeData->neg_elimination_rate_constant, concentration_today, duration);	
    
    // IV factors
    //NOTE: we ignore any potential overlap with remaining concentrations in blood.
    //I.e. if QN is still in blood from previous pills/IV, this factor won't be quite right.
    for( list<IV_dose>::const_iterator dose = IV_doses.begin(); dose != IV_doses.end(); ++dose ){
	//FIXME: calculate new factor
	double factor = 1.0;
	totalFactor *= factor;
    }
    
    return totalFactor;			/* KW -	Returning drug effect per day, per drug */
}

bool LSTMDrug::updateConcentration () {
    double startTime = 0.0;		// as in calculateDrugFactor()
    
    for (multimap<double,double>::const_iterator dose = doses.begin(); dose!=doses.end(); ++dose) {
	if( dose->first >= 1.0 )
	    break;
	
	double duration = dose->first - startTime;
	concentration *= exp(typeData->neg_elimination_rate_constant *  duration);
	concentration += dose->second;
	
	startTime = dose->first;	
    }
   
    double duration = 1.0 - startTime;
    concentration *= exp(typeData->neg_elimination_rate_constant *  duration);
    
    // Clear today's dose list — they've been added to concentration now.
    multimap<double,double>::iterator firstTomorrow = doses.lower_bound( 1.0 );
    doses.erase( doses.begin(), firstTomorrow );
    
    // Now we've removed today's doses, subtract a day from times of tomorrow's doses.
    // Keys are read-only, so we have to create a copy.
    multimap<double,double> newDoses;
    for (multimap<double,double>::const_iterator dose = doses.begin(); dose!=doses.end(); ++dose) {
	// tomorrow's dose; decrease time counter by a day
	newDoses.insert( make_pair<double,double>( dose->first - 1.0, dose->second ) );
    }
    doses.swap( newDoses );	// assign it modified doses (swap may be faster than assign)
    
    // Clear today's IV administrations
    IV_doses.clear();
    
    // return true when concentration is no longer significant:
    return concentration < typeData->negligible_concentration;
}

} }