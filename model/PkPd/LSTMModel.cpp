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

#include "PkPd/LSTMModel.h"
#include "PkPd/LSTMMedicate.h"
#include "PkPd/LSTMTreatments.h"
#include "mon/reporting.h"
#include "util/checkpoint_containers.h"
#include "util/errors.h"

#include "schema/scenario.h"

#include <cassert>

namespace OM { namespace PkPd {

// ———  static init functions  ———

void LSTMModel::init( const scnXml::Scenario& scenario ){
    if (scenario.getPharmacology().present()) {
        LSTMDrugType::init(scenario.getPharmacology().get().getDrugs());
        LSTMTreatments::init(scenario.getPharmacology().get().getTreatments());
    }
}

// ———  non-static set up / tear down functions  ———

void LSTMModel::checkpoint (istream& stream) {
    size_t numDrugs;	// type must be same as m_drugs.size()
    numDrugs & stream;
    util::checkpoint::validateListSize (numDrugs);
    for (size_t i=0; i<numDrugs; ++i) {
        size_t index;
        index & stream;
        m_drugs.push_back (LSTMDrug (LSTMDrugType::getDrug(index)));
        m_drugs.back() & stream;
    }
    medicateQueue & stream;
}

void LSTMModel::checkpoint (ostream& stream) {
    m_drugs.size() & stream;
    for (list<LSTMDrug>::iterator it=m_drugs.begin(); it!=m_drugs.end(); ++it) {
        it->getIndex() & stream;
        (*it) & stream;
    }
    medicateQueue & stream;
}


// ———  non-static simulation time functions  ———

void LSTMModel::prescribe(size_t schedule, size_t dosage, double age, double body_mass){
    const DosageTable& table = dosages[dosage];
    double key = table.useMass ? body_mass : age;
    double doseMult = dosages[dosage].getMultiplier( key );
    foreach( MedicateData& medicateData, schedules[schedule].medications ){
        medicateQueue.push_back( medicateData.multiplied(doseMult) );
    }
}

void LSTMModel::medicate(double body_mass){
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

void LSTMModel::medicateDrug(size_t typeIndex, double qty, double time, double duration, double bodyMass) {
    list<LSTMDrug>::iterator drug = m_drugs.begin();
    while (drug != m_drugs.end()) {
        if (drug->getIndex() == typeIndex)
        goto medicateGotDrug;
        ++drug;
    }
    // No match, so insert one:
    m_drugs.push_front (LSTMDrug(LSTMDrugType::getDrug(typeIndex)));
    drug = m_drugs.begin();	// the drug we just added
    
    medicateGotDrug:
    if( duration > 0.0 ){
        drug->medicateIV (time, duration, qty);
    }else{      // 0 or NaN
        drug->medicate (time, qty, bodyMass);
    }
}

double LSTMModel::getDrugFactor (uint32_t genotype) {
    double factor = 1.0; //no effect
    
    foreach( LSTMDrug& drug, m_drugs ){
        double drugFactor = drug.calculateDrugFactor(genotype);
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
void LSTMModel::decayDrugs () {
    // for each item in m_drugs, remove if DecayPredicate::operator() returns true (so calls decay()):
    m_drugs.remove_if (DecayPredicate());
}

void LSTMModel::summarize(const Host::Human& human) const{
    foreach( const LSTMDrug& drug, m_drugs ){
        assert( drug.getConcentration() > 0 );
        size_t index = drug.getIndex();
        mon::reportMHPI( mon::MHR_HOSTS_POS_DRUG_CONC, human, index, 1 );
        mon::reportMHPF( mon::MHF_LOG_DRUG_CONC, human, index, log(drug.getConcentration()) );
    }
}

} }