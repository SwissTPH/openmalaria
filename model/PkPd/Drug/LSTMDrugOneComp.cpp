/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "PkPd/Drug/LSTMDrugOneComp.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <algorithm>

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrugOneComp::LSTMDrugOneComp(const LSTMDrugType& type) :
    LSTMDrug(),
    typeData(type),
    concentration (0.0)
{}

//LSTMDrugOneComp::~LSTMDrugOneComp(){}

size_t LSTMDrugOneComp::getIndex() const {
    return typeData.getIndex();
}
double LSTMDrugOneComp::getConcentration() const {
    return concentration;
}

void LSTMDrugOneComp::medicate(double time, double qty, double bodyMass)
{
    _medicate(time, qty, typeData.getVolumeOfDistribution() * bodyMass);
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrugOneComp::calculateDrugFactor(uint32_t genotype) {
    /* Survival factor of the parasite (this multiplies the parasite density).
    Calculated below for each time interval. */
    double totalFactor = 1.0;
    
    // Make a copy of concetration and use that over today. Don't adjust concentration because this
    // function may be called multiple times (or not at all) in a day.
    double concentration_today = concentration; // mg / l
    
    const LSTMDrugPD& drugPD = typeData.getPD(genotype);
    
    // Make sure we have a dose at both time 0 and time 1
    //NOTE: this forces function to be non-const and not thread-safe over the same human (probably not an issue)
    //TODO(performance): can we use a faster allocator? Or avoid allocating at all?
    if( doses.begin()->first != 0.0 ){
        doses.insert( doses.begin(), make_pair( 0.0, DoseParams() ) );
    }
    if( doses.count( 1.0 ) == 0 ){
        doses.insert( make_pair( 1.0, DoseParams() ) );
    }
    
    DoseMap::const_iterator dose = doses.begin();
    DoseMap::const_iterator next_dose = dose;
    ++next_dose;
    while (next_dose!=doses.end()) {
        double time_to_next = next_dose->first - dose->first;
        if( dose->second.duration == 0.0 ){
            // Oral dose
            concentration_today += dose->second.qty;
            
            totalFactor *= drugPD.calcFactor( typeData, concentration_today, time_to_next );
        } else {
            // IV dose
            assert( util::vectors::approxEqual(time_to_next, dose->second.duration) );
            
            totalFactor *= drugPD.calcFactorIV( typeData, concentration_today, time_to_next, dose->second.qty );
        }
        
        dose = next_dose;
        ++next_dose;
        if( dose->first >= 1.0 )
            break;      // we know this and any more doses happen tomorrow; don't calculate factors now
    }
    
    return totalFactor; // Drug effect per day per drug per parasite
}

bool LSTMDrugOneComp::updateConcentration () {
    // Make sure we have a dose at both time 0 and time 1
    //TODO(performance): can we use a faster allocator? Or avoid allocating at all?
    if( doses.begin()->first != 0.0 ){
        doses.insert( doses.begin(), make_pair( 0.0, DoseParams() ) );
    }
    if( doses.count( 1.0 ) == 0 ){
        doses.insert( make_pair( 1.0, DoseParams() ) );
    }
    
    DoseMap::const_iterator dose = doses.begin();
    DoseMap::const_iterator next_dose = dose;
    ++next_dose;
    while (next_dose!=doses.end()) {
        double time_to_next = next_dose->first - dose->first;
        if( dose->second.duration == 0.0 ){
            // Oral dose
            concentration += dose->second.qty;
            typeData.updateConcentration( concentration, time_to_next );
        } else {
            // IV dose
            assert( util::vectors::approxEqual(time_to_next, dose->second.duration) );
            
            typeData.updateConcentrationIV( concentration, time_to_next, dose->second.qty );
        }
        
        dose = next_dose;
        ++next_dose;
        if( dose->first >= 1.0 )
            break;      // we know this and any more doses happen tomorrow; don't calculate factors now
    }
    
    // Clear today's dose list — they've been added to concentration now.
    DoseMap::iterator firstTomorrow = doses.lower_bound( 1.0 );
    doses.erase( doses.begin(), firstTomorrow );
    
    //TODO(performance): is there some way we can avoid copying here? Cache all possible simulated days?
    // Now we've removed today's doses, subtract a day from times of tomorrow's doses.
    // Keys are read-only, so we have to create a copy.
    DoseMap newDoses;
    for (DoseMap::const_iterator dose = doses.begin(); dose!=doses.end(); ++dose) {
        // tomorrow's dose; decrease time counter by a day
        newDoses.insert( make_pair( dose->first - 1.0, dose->second ) );
    }
    doses.swap( newDoses );     // assign it modified doses (swap may be faster than assign)
    
    util::streamValidate( concentration );
    
    // return true when concentration is no longer significant:
    return concentration < typeData.getNegligibleConcentration();
}

void LSTMDrugOneComp::checkpoint(ostream& stream){
    concentration & stream;
}
void LSTMDrugOneComp::checkpoint(istream& stream){
    concentration & stream;
}

}
}