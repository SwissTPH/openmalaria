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

#include "PkPd/HoshenPkPdModel.h"

namespace OM { namespace PkPd {
    
// -----  static functions  -----

void HoshenPkPdModel::init() {
}


// -----  non-static set up / tear down functions  -----

HoshenPkPdModel::HoshenPkPdModel () {}
HoshenPkPdModel::~HoshenPkPdModel () {}

void HoshenPkPdModel::checkpoint (istream& stream) {
    size_t numDrugs;	// type must be same as _drugs.size()
    numDrugs & stream;
    validateListSize (numDrugs);
    for (size_t i=0; i<numDrugs; i++) {
	string abbrev;
	abbrev & stream;
	_drugs.push_back (HoshenDrug ((HoshenDrugType*)DrugType::getDrug(abbrev)));
	_drugs.back() & stream;
    }
}

void HoshenPkPdModel::checkpoint (ostream& stream) {
    _drugs.size() & stream;
    for (list<HoshenDrug>::iterator it=_drugs.begin(); it!=_drugs.end(); it++) {
	it->getAbbreviation() & stream;
	(*it) & stream;
    }
}


// -----  non-static simulation time functions  -----

void HoshenPkPdModel::medicate(string drugAbbrev, double qty, int time, double ageYears) {
  list<HoshenDrug>::iterator drug = _drugs.begin();
  while (drug != _drugs.end()) {
    if (drug->getAbbreviation() == drugAbbrev)
      goto medicateGotDrug;
    ++drug;
  }
  // No match, so insert one:
  _drugs.push_front (HoshenDrug((HoshenDrugType*)DrugType::getDrug(drugAbbrev)));
  drug = _drugs.begin();	// the drug we just added
  
  medicateGotDrug:
  double weight = ageToWeight(ageYears);
  drug->addDose (qty*drug->getAbsorptionFactor()/weight, time);
}

struct DecayPredicate {
  bool operator() (HoshenDrug& drug) {
    return drug.decay ();
  }
};
void HoshenPkPdModel::decayDrugs (double) {
  // for each item in _drugs, remove if DecayPredicate::operator() returns true (so calls decay()):
  _drugs.remove_if (DecayPredicate());
}

double HoshenPkPdModel::getDrugFactor (uint32_t proteome_ID, double ageYears) {
  // We will choose for now the smallest factor (ie, most impact)
  double factor = 1.0; //no effect
  
  for (list<HoshenDrug>::const_iterator it=_drugs.begin(); it!=_drugs.end(); it++) {
    double drugFactor = it->calculateDrugFactor(proteome_ID);
    if (drugFactor < factor) {
      factor = drugFactor;
    }
  }
  return factor;
}

} }