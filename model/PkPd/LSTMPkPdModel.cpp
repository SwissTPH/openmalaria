/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "PkPd/LSTMPkPdModel.h"

namespace OM { namespace PkPd {
    
// -----  static functions  -----

void LSTMPkPdModel::init() {
}


// -----  non-static set up / tear down functions  -----

LSTMPkPdModel::LSTMPkPdModel () {
    // TODO (LSTM): can add initialization, heterogeneity, etc., here
    //metabolismMultiplier = 3;
}
LSTMPkPdModel::~LSTMPkPdModel () {}

void LSTMPkPdModel::checkpoint (istream& stream) {
    size_t numDrugs;	// type must be same as _drugs.size()
    numDrugs & stream;
    validateListSize (numDrugs);
    for (size_t i=0; i<numDrugs; i++) {
	string abbrev;
	abbrev & stream;
	_drugs.push_back (LSTMDrug (LSTMDrugType::getDrug(abbrev)));
	_drugs.back() & stream;
    }
}

void LSTMPkPdModel::checkpoint (ostream& stream) {
    _drugs.size() & stream;
    for (list<LSTMDrug>::iterator it=_drugs.begin(); it!=_drugs.end(); it++) {
	it->getAbbreviation() & stream;
	(*it) & stream;
    }
}


// -----  non-static simulation time functions  -----

void LSTMPkPdModel::medicate(string drugAbbrev, double qty, double time, double ageYears) {
  list<LSTMDrug>::iterator drug = _drugs.begin();
  while (drug != _drugs.end()) {
    if (drug->getAbbreviation() == drugAbbrev)
      goto medicateGotDrug;
    ++drug;
  }
  // No match, so insert one:
  _drugs.push_front (LSTMDrug(LSTMDrugType::getDrug(drugAbbrev)));
  drug = _drugs.begin();	// the drug we just added
  
  medicateGotDrug:
  drug->medicate (time, qty, ageToWeight (ageYears));
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

double LSTMPkPdModel::getDrugFactor (uint32_t proteome_ID) {
    double factor = 1.0; //no effect
    
    for (list<LSTMDrug>::iterator it=_drugs.begin(); it!=_drugs.end(); it++) {
	double drugFactor = it->calculateDrugFactor(proteome_ID);
	factor *= drugFactor;
    }
    return factor;
}

uint32_t LSTMPkPdModel::new_proteome_ID () {
    return LSTMDrugType::new_proteome_ID ();
}

} }