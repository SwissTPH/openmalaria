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

#include "PkPd/Drug/DrugType.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

/*
 * Static variables and functions
 */

map<string,DrugType*> DrugType::available; 



void DrugType::init () {
}


void DrugType::addDrug(DrugType* drug) {
  string abbrev = drug->abbreviation;
  // Check drug doesn't already exist
    if (available.find (abbrev) != available.end())
    throw invalid_argument (string ("Drug already in registry: ").append(abbrev));
  
  available.insert (pair<string,DrugType*>(abbrev, drug));
}

const DrugType* DrugType::getDrug(string _abbreviation) {
  map<string,DrugType*>::const_iterator i = available.find (_abbreviation);
  if (i == available.end())
    throw OM::util::errors::xml_scenario_error (string ("prescribed non-existant drug ").append(_abbreviation));
  
  return i->second;
}


// -----  Non-static DrugType functions  -----

DrugType::DrugType (string _name, string _abbreviation) {
  name = _name;
  abbreviation = _abbreviation;
}
DrugType::~DrugType () {}

/* Checkpointing functions, which we shouldn't need now. If they are needed:
DrugType::DrugType (istream& in) {
  int num;
  in >> abbreviation;
  in >> name;
  in >> absorptionFactor;
  in >> halfLife; 
  in >> num;
  Global::validateListSize (num);
  for (int i=0; i<num; i++) {
    vector <Mutation*> rule = vector<Mutation*>();
    int numMuts;
    in >> numMuts;
    Global::validateListSize (numMuts);
    for (int j=0; j<numMuts; j++) {
      string proteinName;
      int position;
      char allele;
      in >> proteinName;
      in >> position;
      in >> allele;
      rule.push_back(
          ProteomeManager::getMutation(proteinName, position, allele));
    }
    requiredMutations.push_back(rule);
  }
  in >> num;
  Global::validateListSize (num);
  for (int i=0; i<num; i++) {
    double parameter;
    in >> parameter;
    pdParameters.push_back(parameter);
  }
  in >> num;
  Global::validateListSize (num);
  for (int i=0; i<num; i++) {
    int proteomeID;
    double parameter;
    in >> proteomeID;
    in >> parameter;
    proteomePDParameters[proteomeID] = parameter;
  }
}
void DrugType::write (ostream& out) const {
  out << abbreviation << endl;
  out << name << endl;
  out << absorptionFactor << endl;
  out << halfLife << endl; 
  out << requiredMutations.size() << endl;
  for (vector< vector<Mutation*> >::const_iterator itRM=requiredMutations.begin(); itRM!=requiredMutations.end(); itRM++) {
    vector<Mutation*>::const_iterator itM;
    out << (*itRM).size() << endl;
    for (itM=(*itRM).begin(); itM!=(*itRM).end(); itM++) {
      out << (*itM)->getProteinName() << endl;
      out << (*itM)->getPosition() << endl;
      out << (*itM)->getAllele() << endl;
    }
  }
  out << pdParameters.size() << endl;
  for (vector<double>::const_iterator itPP=pdParameters.begin(); itPP!=pdParameters.end(); itPP++) {
    out << *itPP << endl ;
  }
  out << proteomePDParameters.size() << endl;
  for (map<int,double>::const_iterator itPPD=proteomePDParameters.begin(); itPPD!=proteomePDParameters.end(); itPPD++) {
    out << (*itPPD).first << endl ;
    out << (*itPPD).second << endl ;
  }
}
*/



