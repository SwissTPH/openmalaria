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

#include "PkPd/Drug/LSTMDrugType.h"

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

void LSTMDrugType::init () {
  DrugType::init();
  //TODO: add drugs from scenario file
  LSTMDrugType* s;
  s = new LSTMDrugType("Mefloquine", "MF");
  DrugType::addDrug(s);
  //s = new LSTMDrugType("Chloroquine", "CQ"); - This will be moved to XML later
  //DrugType::addDrug(s);
}

// -----  Non-static DrugType functions  -----

LSTMDrugType::LSTMDrugType (string _name, string _abbreviation)
: DrugType(_name, _abbreviation)
{
	//TODO: temporarily hard-coding values; need to get from XML
	if(_abbreviation == "MF") {
		parameters.max_killing_rate = 3.45;
		parameters.IC50 = 0.6654;
		parameters.slope = 2.5;
		parameters.elimination_rate_constant = 0.036;
		parameters.vol_dist = 20.8;
	} //else if(_abbreviation == "CQ"){...} - This will be moved to XML later
}
LSTMDrugType::~LSTMDrugType () {}


void LSTMDrugType::addPDRule(vector<Mutation*> ruleRequiredMutations, double pdFactor) {
    /*FIXME (possibly totally unwanted)
  requiredMutations.push_back(ruleRequiredMutations);
  pdParameters.push_back(pdFactor);*/
}

void LSTMDrugType::parseProteomeInstances() {
    /*FIXME (possibly totally unwanted)
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
  */
}

} }