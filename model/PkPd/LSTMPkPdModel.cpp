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

#include "PkPd/LSTMPkPdModel.h"
#include "PkPd/LSTMMedicate.h"
#include "util/checkpoint_containers.h"
#include "util/errors.h"

#include "schema/scenario.h"

#include <cassert>

namespace OM { namespace PkPd {

// ———  non-static set up / tear down functions  ———

LSTMPkPdModel::LSTMPkPdModel () {}
LSTMPkPdModel::~LSTMPkPdModel () {}

void LSTMPkPdModel::checkpoint (istream& stream) {
    size_t numDrugs;	// type must be same as _drugs.size()
    numDrugs & stream;
    util::checkpoint::validateListSize (numDrugs);
    for (size_t i=0; i<numDrugs; ++i) {
        size_t index;
        index & stream;
        _drugs.push_back (LSTMDrug (LSTMDrugType::getDrug(index)));
        _drugs.back() & stream;
    }
    medicateQueue & stream;
}

void LSTMPkPdModel::checkpoint (ostream& stream) {
    _drugs.size() & stream;
    for (list<LSTMDrug>::iterator it=_drugs.begin(); it!=_drugs.end(); ++it) {
        it->getIndex() & stream;
        (*it) & stream;
    }
    medicateQueue & stream;
}


// ———  non-static simulation time functions  ———

void LSTMPkPdModel::prescribe(size_t schedule, size_t dosage, double age, double body_mass){
    const DosageTable& table = dosages[dosage];
    double key = table.useMass ? body_mass : age;
    double doseMult = dosages[dosage].getMultiplier( key );
    foreach( MedicateData& medicateData, schedules[schedule].medications ){
        medicateQueue.push_back( medicateData.multiplied(doseMult) );
    }
}

void LSTMPkPdModel::medicate(double body_mass){
    if( medicateQueue.empty() ) return;
    
    // Process pending medications (in interal queue) and apply/update:
    list<MedicateData>::iterator it = medicateQueue.begin();
    while( it != medicateQueue.end() ){
        if( it->time < 1.0 ){   // Medicate medications to be prescribed starting at the next time-step
            // This function could be inlined, except for uses in testing:
            medicateDrug (it->drug, it->qty, it->time, it->duration, body_mass);
            it = medicateQueue.erase (it);
        }else{  // and decrement treatment seeking delay for the rest
            it->time -= 1.0;
            ++it;
        }
    }
}

void LSTMPkPdModel::medicateDrug(size_t typeIndex, double qty, double time, double duration, double bodyMass) {
    list<LSTMDrug>::iterator drug = _drugs.begin();
    while (drug != _drugs.end()) {
        if (drug->getIndex() == typeIndex)
        goto medicateGotDrug;
        ++drug;
    }
    // No match, so insert one:
    _drugs.push_front (LSTMDrug(LSTMDrugType::getDrug(typeIndex)));
    drug = _drugs.begin();	// the drug we just added
    
    medicateGotDrug:
    if( duration > 0.0 ){
        drug->medicateIV (time, duration, qty);
    }else{      // 0 or NaN
        drug->medicate (time, qty, bodyMass);
    }
}

double LSTMPkPdModel::getDrugFactor (uint32_t genotype) {
    double factor = 1.0; //no effect
    
    for (list<LSTMDrug>::iterator it=_drugs.begin(); it!=_drugs.end(); ++it) {
        double drugFactor = it->calculateDrugFactor(genotype);
        factor *= drugFactor;
    }
    return factor;
}

// This may look complicated but its just some machinery to call updateConcentration() and return its result
class DecayPredicate {
public:
    bool operator() (LSTMDrug& drug) {
        return drug.updateConcentration();
    }
};
void LSTMPkPdModel::decayDrugs () {
    // for each item in _drugs, remove if DecayPredicate::operator() returns true (so calls decay()):
    _drugs.remove_if (DecayPredicate());
}

} }