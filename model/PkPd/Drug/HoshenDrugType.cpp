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

#include "PkPd/Drug/HoshenDrugType.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

namespace OM { namespace PkPd {
    
/*
 * Static variables and functions
 */

//map<string,HoshenDrugType*> HoshenDrugType::available; 

void HoshenDrugType::init () {
  DrugType::init();
  Mutation* crt76 = ProteomeManager::getMutation(string("CRT"), 76, 'T');
  HoshenDrugType* s;
  //s = new DrugType("Sulfadoxine", "S", 0.1, 10*24*60); //Invented values
  //DrugType::addDrug(s);
  s = new HoshenDrugType("Chloroquine", "CQ", 0.02, 45*24*60); //Based on Hoshen
  vector<Mutation*> crt76L;
  crt76L.push_back(crt76);
  s->addPDRule(crt76L, 204.0);
  s->addPDRule(vector<Mutation*>(), 68.0);
  s->parseProteomeInstances();
  DrugType::addDrug(s);
}

// -----  Non-static DrugType functions  -----

HoshenDrugType::HoshenDrugType (string _name, string _abbreviation,
    double _absorptionFactor, double _halfLife)
: DrugType(_abbreviation)
{
  //name = _name;
  //abbreviation = _abbreviation;
  absorptionFactor = _absorptionFactor;
  halfLife = _halfLife;
}
HoshenDrugType::~HoshenDrugType () {}


void HoshenDrugType::addPDRule(vector<Mutation*> ruleRequiredMutations, double pdFactor) {
  requiredMutations.push_back(ruleRequiredMutations);
  pdParameters.push_back(pdFactor);
}

void HoshenDrugType::parseProteomeInstances() {
  vector<ProteomeInstance> instances = ProteomeInstance::getInstances();
  int numRules = requiredMutations.size();
  for (vector<ProteomeInstance>::const_iterator it=instances.begin(); it !=instances.end(); it++) {
    //cerr << " Here goes instance";
    for(int rule=0; rule<numRules; rule++) {
      if (it->hasMutations(requiredMutations[rule])) {
	proteomePDParameters[it->getProteomeID()] = pdParameters[rule];
	//cerr << " rule: " << rule << "\n";
	break;
      }
    }
  }
}

} }