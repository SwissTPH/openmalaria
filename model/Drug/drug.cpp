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

#include "Drug/drug.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

/*
 * Static variables and functions
 */

map<string,DrugType> DrugType::available; 
double Drug::minutesPerTimeStep;


void DrugType::init () {
  Drug::init ();
  Mutation* crt76 = ProteomeManager::getMutation(string("CRT"), 76, 'T');
  DrugType* s;
  //s = new DrugType("Sulfadoxine", "S", 0.1, 10*24*60); //Invented values
  //DrugType::addDrug(s);
  s = new DrugType("Chloroquine", "CQ", 0.02, 45*24*60); //Based on Hoshen
  vector<Mutation*> crt76L;
  crt76L.push_back(crt76);
  s->addPDRule(crt76L, 204.0);
  s->addPDRule(vector<Mutation*>(), 68.0);
  s->parseProteomeInstances();
  DrugType::addDrug(s);
}
void Drug::init () {
  minutesPerTimeStep = Global::interval * 24*60;
}

void DrugType::addDrug(DrugType* drug) {
  string abbrev = drug->abbreviation;
  // Check drug doesn't already exist
  if (available.find (abbrev) != available.end())
    throw invalid_argument (string ("Drug already in registry: ").append(abbrev));
  
  available.insert (pair<string,DrugType>(abbrev, *drug));
}

const DrugType* DrugType::getDrug(string _abbreviation) {
  map<string,DrugType>::const_iterator i = available.find (_abbreviation);
  if (i == available.end())
    throw xml_scenario_error (string ("prescribed non-existant drug ").append(_abbreviation));
  
  return &i->second;
}


// -----  Non-static DrugType functions  -----

DrugType::DrugType (string _name, string _abbreviation,
    double _absorptionFactor, double _halfLife)
{
  name = _name;
  abbreviation = _abbreviation;
  absorptionFactor = _absorptionFactor;
  halfLife = _halfLife;
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


void DrugType::addPDRule(vector<Mutation*> ruleRequiredMutations, double pdFactor) {
  requiredMutations.push_back(ruleRequiredMutations);
  pdParameters.push_back(pdFactor);
}

void DrugType::parseProteomeInstances() {
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


// -----  Dose methods  -----

Dose::Dose (istream& in) {
  in >> x;
  in >> y;
}
void Dose::write (ostream& out) const {
  out << x << endl;
  out << y << endl;
}

/*
 * Drug code
 */

Drug::Drug(const DrugType* type) :
  typeData (type),
  _concentration (0),
  _nextConcentration (0)
{}
Drug::Drug (const DrugType* type, istream& in) :
  typeData (type)
{
  in >> _concentration;
  // No need to checkpoint; just recalculate _nextConcentration
  _nextConcentration = _concentration * decayFactor (minutesPerTimeStep);
  int num;
  in >> num;
  Global::validateListSize (num);
  for (int i = 0; i < num; ++i)
    doses.push_back (Dose (in));
}

void Drug::write (ostream& out) const {
  out << _concentration << endl;
  assert (doses.size() == 0);	// doses not used yet (remove this eventually)
  out << doses.size() << endl;
  for (deque<Dose>::const_iterator it = doses.begin(); it != doses.end(); ++it)
    it->write (out);
}


void Drug::addDose (double concentration, int delay) {
  if (delay == 0) {
    _concentration += concentration;
    _nextConcentration = _concentration * decayFactor (minutesPerTimeStep);
  } else {
    assert (false);	//NOTE: Not properly dealt with yet (as it is only relevant for ACTs).
    // Only adding doses for this timestep is supported
    assert (delay>0 && delay<minutesPerTimeStep);
    double nextConcentration = concentration*decayFactor (minutesPerTimeStep-delay);
    doses.push_back (Dose (nextConcentration, 0 /*FIXME*/));
  }
}

double Drug::calculateDrugFactor(ProteomeInstance* infProteome) const {
  //Returning an average of 2 points
  double param = typeData->proteomePDParameters.find(infProteome->getProteomeID())->second;
  double startFactor = 3.8/(1+param/_concentration);
  double endFactor = 3.8/(1+param/_nextConcentration);
  return exp(-(startFactor + endFactor)/2);
}

double Drug::decayFactor (double time) {
  //k = log(2)/halfLife
  return exp(-time*log(2.0)/typeData->halfLife);
}

bool Drug::decay() {
  _concentration = _nextConcentration;
  _nextConcentration = _concentration * decayFactor (minutesPerTimeStep);
  //TODO: if concentration is negligible, return true to clean up this object
  return false;
}
