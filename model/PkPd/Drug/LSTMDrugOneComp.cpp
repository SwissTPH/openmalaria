/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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
double LSTMDrugOneComp::calculateDrugFactor(uint32_t genotype) const {
    /* Survival factor of the parasite (this multiplies the parasite density).
    Calculated below for each time interval. */
    double totalFactor = 1.0;
    
    // Make a copy of concetration and use that over today. Don't adjust concentration because this
    // function may be called multiple times (or not at all) in a day.
    double concentration_today = concentration; // mg / l
    
    const LSTMDrugPD& drugPD = typeData.getPD(genotype);
    
    double time = 0.0;
    typedef pair<double,double> TimeConc;
    foreach( TimeConc time_conc, doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( time < time_conc.first ){
                totalFactor *= drugPD.calcFactor( typeData, concentration_today, time_conc.first - time );
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add dose (instantaneous absorbtion):
            concentration_today += time_conc.second;
        }else/*i.e. tomorrow or later*/{
            break;
        }
    }
    if( time < 1.0 ){
        totalFactor *= drugPD.calcFactor( typeData, concentration_today, 1.0 - time );
    }
    
    return totalFactor; // Drug effect per day per drug per parasite
}

bool LSTMDrugOneComp::updateConcentration () {
    double time = 0.0;
    size_t doses_taken = 0;
    typedef pair<double,double> TimeConc;
    foreach( TimeConc& time_conc, doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( time < time_conc.first ){
                // exponential decay of drug concentration:
                concentration *= exp(typeData.getNegElimintationRateConst() * (time_conc.first - time));
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add dose (instantaneous absorbtion):
            concentration += time_conc.second;
            doses_taken += 1;
        }else /*i.e. tomorrow or later*/{
            time_conc.first -= 1.0;
        }
    }
    if( time < 1.0 ){
        // exponential decay of drug concentration:
        concentration *= exp(typeData.getNegElimintationRateConst() * (1.0 - time));
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    doses.erase(doses.begin(), doses.begin() + doses_taken);
    
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