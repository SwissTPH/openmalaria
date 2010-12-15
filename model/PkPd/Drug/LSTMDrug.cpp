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
#include "util/StreamValidator.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

#include <gsl/gsl_integration.h>

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrug::LSTMDrug(const LSTMDrugType& type) :
    typeData (&type),
    concentration (0.0)
{}


// Create our list of doses. Optimise for the case where we only have 1 or less per day, but be able to handle more.
// If two oral doses coincide, they can be combined but doing so is not essential.
// End of an IV dose can be represented by an oral dose of 0 qty potentially.
// IV doses must be split over the end of the day and if an oral dose occurs in the middle.
// Overlapping IV doses are not supported.

void LSTMDrug::medicate (double time, double qty, double weight) {
    double conc = qty / (typeData->vol_dist * weight);
    // multimap insertion: is ordered
    DoseMap::iterator lastInserted =
    doses.insert (doses.end(), make_pair (time, DoseParams( conc, 0 )));
    check_split_IV( lastInserted );
}

void LSTMDrug::medicateIV (double duration, double endTime, double qty ) {
    assert( duration > 0.0 );
    
    double infusRate = qty / duration;	// mg/day
    DoseMap::iterator lastInserted =
    doses.insert(doses.end(), make_pair( endTime-duration, DoseParams( infusRate, duration ) ) );    
    check_split_IV( lastInserted );
    
    doses.insert( doses.end(), make_pair( endTime, DoseParams() ) );
}

void LSTMDrug::check_split_IV( DoseMap::iterator lastInserted ){
    for( DoseMap::iterator it = doses.begin(); it != doses.end(); ++it ){
        if( it == lastInserted )
            continue;
        
        if( it->first == lastInserted->first ){
            if( it->second.duration == 0.0 && lastInserted->second.duration == 0.0 ){
                it->second.qty += lastInserted->second.qty;
                doses.erase( lastInserted );
                return;
            } else if( it->second.duration != 0.0 && lastInserted->second.duration != 0.0 ){
                throw util::xml_scenario_error( "IV/oral medications overlap — not supported!" );
            }
        } else if( it->first < lastInserted->first ){
            if( it->first + it->second.duration > lastInserted->first ){
                // TODO: if one is an oral dose, could split IV
                throw util::xml_scenario_error( "IV/oral medications overlap — not supported!" );
            }
        } else if( it->first > lastInserted->first ){
            if( lastInserted->first + lastInserted->second.duration > it->first ){
                // TODO: if one is an oral dose, could split IV
                throw util::xml_scenario_error( "IV/oral medications overlap — not supported!" );
            }
        }
    }
    
    // Make sure one dose is evaluated separately on each day
    double dayEnd = 1.0;
    while( lastInserted->first + lastInserted->second.duration > dayEnd ){
        double remaining_duration = lastInserted->first + lastInserted->second.duration - dayEnd;
        lastInserted->second.duration = dayEnd - lastInserted->first;
        lastInserted = doses.insert( doses.end(), make_pair( dayEnd, DoseParams( lastInserted->second.qty, remaining_duration ) ) );
        dayEnd += 1.0;
    }
}


struct IV_conc_params {
    double C0;
    double ivRate; // dose->second.qty
    double neg_elimination_rate_constant;
    double elim_rate_dist;      // -neg_elimination_rate_constant * vol_dist
    double slope;
    double max_kill_rate;
    double IC50_pow_slope;
};
/** Function for calculating killing factor of IV and requiring integration.
 */
double func_IV_conc( double t, void* pp ){
    IV_conc_params *p = static_cast<IV_conc_params*>( pp );
    
    double conc_decay = exp(p->neg_elimination_rate_constant * t);
    double infusion = p->ivRate * (1.0 - conc_decay ) / p->elim_rate_dist
        + p->C0 * conc_decay;
    double infusion_pow_slope = pow(infusion, p->slope);
    double fC = p->max_kill_rate * infusion_pow_slope / infusion_pow_slope + p->IC50_pow_slope;
    return fC;
}


// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrug::calculateDrugFactor(uint32_t proteome_ID) {
    /* Survival factor of the parasite (this multiplies the parasite density).
    Calculated below for each time interval. */
    double totalFactor = 1.0;
    
    // Make a copy of concetration and use that over today. Don't adjust concentration because this
    // function may be called multiple times (or not at all) in a day.
    double concentration_today = concentration;
    
    size_t allele = (proteome_ID >> typeData->allele_rshift) & typeData->allele_mask;
    const LSTMDrugPDParameters& PD_params = typeData->PD_params[allele];
    
    // Make sure we have a dose at both time 0 and time 1
    //TODO: analyse efficiency of this method
    //NOTE: this forces function to be non-const and not thread-safe over the same human (probably not an issue)
    if( doses.begin()->first != 0.0 ){
        doses.insert( doses.begin(), make_pair( 0.0, DoseParams() ) );
    }
    if( doses.count( 1.0 ) == 0 ){
        doses.insert( make_pair( 1.0, DoseParams() ) );
    }
    
    DoseMap::const_iterator dose = doses.begin();
    DoseMap::const_iterator next_dose = dose;
    next_dose++;
    while (next_dose!=doses.end()) {
        double duration = next_dose->first - dose->first;
	if( dose->second.duration == 0.0 ){
            // Oral dose
            concentration_today += dose->second.qty;
            
            double initial_conc = concentration_today;
            concentration_today *= exp(typeData->neg_elimination_rate_constant *  duration);
            
            // Note: these look a little different from original equations because PD_params.IC50_pow_slope
            // and PD_params.power are calculated when read from the scenario document instead of here.
            double numerator = PD_params.IC50_pow_slope + pow(concentration_today,PD_params.slope);
            double denominator = PD_params.IC50_pow_slope + pow(initial_conc,PD_params.slope);
            totalFactor *= pow( numerator / denominator, PD_params.power );
        } else {
            // IV dose
            assert( duration == dose->second.duration );
            
            //TODO: implement caching here: if already evaluated for this dose
            // and these PD_params, then reuse result.
            
            IV_conc_params p;
            p.C0 = concentration_today;
            p.ivRate =  dose->second.qty;
            p.neg_elimination_rate_constant = typeData->neg_elimination_rate_constant;
            p.elim_rate_dist =-p.neg_elimination_rate_constant * typeData->vol_dist;
            p.slope = PD_params.slope;
            p.max_kill_rate = PD_params.max_killing_rate;
            p.IC50_pow_slope = PD_params.IC50_pow_slope;
            
            gsl_function F;
            F.function = &func_IV_conc;
            F.params = static_cast<void*>(&p);
            
            //TODO: consider desired limits
            double abs_eps = 0.0, rel_eps = 1e-7;
            double intfC, err_eps;
            size_t n_evals;
            
            gsl_integration_qng (&F, 0.0, duration, abs_eps, rel_eps, &intfC, &err_eps, &n_evals); 
            
#ifndef NDEBUG
            cout << "IV: iterations: "<<n_evals<<"; error epsilon: "<<err_eps<<endl;
#endif
            //TODO: check err_eps is OK
            //TODO: check qng is the best integration algorithm
            
            totalFactor *= 1.0 / exp( intfC );
            
            concentration_today *= exp(typeData->neg_elimination_rate_constant * duration);
            concentration_today += dose->second.qty
                * (1.0 - exp(typeData->neg_elimination_rate_constant * duration) )
                / p.elim_rate_dist;
        }
        
        dose = next_dose;
        next_dose++;
        if( dose->first >= 1.0 )
            break;      // we know this and any more doses happen tomorrow; don't calculate factors now
    }
    
    return totalFactor;	// Drug effect per day per drug per parasite
}

bool LSTMDrug::updateConcentration () {
    // Make sure we have a dose at both time 0 and time 1
    //TODO: analyse efficiency of this method
    if( doses.begin()->first != 0.0 ){
        doses.insert( doses.begin(), make_pair( 0.0, DoseParams() ) );
    }
    if( doses.count( 1.0 ) == 0 ){
        doses.insert( make_pair( 1.0, DoseParams() ) );
    }
    
    DoseMap::const_iterator dose = doses.begin();
    DoseMap::const_iterator next_dose = dose;
    next_dose++;
    while (next_dose!=doses.end()) {
        double duration = next_dose->first - dose->first;
        if( dose->second.duration == 0.0 ){
            // Oral dose
            concentration += dose->second.qty;
            concentration *= exp(typeData->neg_elimination_rate_constant *  duration);
        } else {
            // IV dose
            assert( duration == dose->second.duration );
            
            //NOTE: could calculate equivalent oral dose, add at start, then decay everything instead
            // then decay could be done irrespective of case
            concentration *= exp(typeData->neg_elimination_rate_constant *  duration);
            concentration += dose->second.qty
                * (1.0 - exp(typeData->neg_elimination_rate_constant * duration) )
                / ( -typeData->neg_elimination_rate_constant * typeData->vol_dist );
        }
        
        dose = next_dose;
        next_dose++;
        if( dose->first >= 1.0 )
            break;      // we know this and any more doses happen tomorrow; don't calculate factors now
    }
    
    // Clear today's dose list — they've been added to concentration now.
    DoseMap::iterator firstTomorrow = doses.lower_bound( 1.0 );
    doses.erase( doses.begin(), firstTomorrow );
    
    // Now we've removed today's doses, subtract a day from times of tomorrow's doses.
    // Keys are read-only, so we have to create a copy.
    DoseMap newDoses;
    for (DoseMap::const_iterator dose = doses.begin(); dose!=doses.end(); ++dose) {
	// tomorrow's dose; decrease time counter by a day
	newDoses.insert( make_pair<double,DoseParams>( dose->first - 1.0, dose->second ) );
    }
    doses.swap( newDoses );	// assign it modified doses (swap may be faster than assign)
    
    util::streamValidate( concentration );
    
    // return true when concentration is no longer significant:
    return concentration < typeData->negligible_concentration;
}

}

namespace util { namespace checkpoint {
    
void operator& (multimap<double,PkPd::DoseParams> x, ostream& stream) {
    x.size() & stream;
    for (multimap<double,PkPd::DoseParams>::iterator pos = x.begin (); pos != x.end() ; ++pos) {
        pos->first & stream;
        pos->second & stream;
    }
}
void operator& (multimap<double,PkPd::DoseParams>& x, istream& stream) {
    size_t l;
    l & stream;
    validateListSize (l);
    x.clear ();
    multimap<double,PkPd::DoseParams>::iterator pos = x.begin ();
    for (size_t i = 0; i < l; ++i) {
        double s;
        PkPd::DoseParams t;
        s & stream;
        t & stream;
        pos = x.insert (pos, make_pair (s,t));
    }
    assert( x.size() == l );
}

} }

}