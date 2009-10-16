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

#include "Drug/PkPdDrug.h"
#include "Drug/drug.h"

// -----  static functions  -----

void PkPdDrug::init() {
  DrugType::init();
}


// -----  non-static set up / tear down functions  -----

PkPdDrug::PkPdDrug () {}
PkPdDrug::~PkPdDrug () {}

PkPdDrug::PkPdDrug (istream& in) {
  int numDrugs;
  in >> numDrugs;
  Global::validateListSize (numDrugs);
  for (int i=0; i<numDrugs; i++) {
    string abbrev;
    in >> abbrev;
    _drugs.push_back (Drug (DrugType::getDrug(abbrev), in));
  }
}
void PkPdDrug::write (ostream& out) {
  out << _drugs.size() << endl;
  for (list<Drug>::const_iterator it=_drugs.begin(); it!=_drugs.end(); it++) {
    out << it->getAbbreviation() << endl;
    it->write (out);
  }
}


// -----  non-static simulation time functions  -----

void PkPdDrug::medicate(string drugAbbrev, double qty, int time, double age, double weight) {
  list<Drug>::iterator drug = _drugs.begin();
  while (drug != _drugs.end()) {
    if (drug->getAbbreviation() == drugAbbrev)
      goto medicateGotDrug;
    ++drug;
  }
  // No match, so insert one:
  _drugs.push_front (Drug(DrugType::getDrug(drugAbbrev)));
  drug = _drugs.begin();	// the drug we just added
  
  medicateGotDrug:
  drug->addDose (qty*drug->getAbsorptionFactor()/weight, time);
}

struct DecayPredicate {
  bool operator() (Drug& drug) {
    return drug.decay ();
  }
};
void PkPdDrug::decayDrugs () {
  // for each item in _drugs, remove if DecayPredicate::operator() returns true (so calls decay()):
  _drugs.remove_if (DecayPredicate());
}

double PkPdDrug::getDrugFactor (const ProteomeInstance* infProteome) {
  // We will choose for now the smallest (ie, most impact)
  
  double factor = 1.0; //no effect
  for (list<Drug>::const_iterator it=_drugs.begin(); it!=_drugs.end(); it++) {
    double drugFactor = it->calculateDrugFactor(infProteome);
    if (drugFactor < factor) {
      factor = drugFactor;
    }
  }
  return factor;
}
