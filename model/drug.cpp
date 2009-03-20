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

#include <iostream>
#include <cmath>
#include "Infection.h"
#include "human.h"
#include "drug.h"

using namespace std;

/*
 * Private variables
 */

static int withinHostTimestep;
static int simulatorTimestep;
static ProteomeManager* manager;

/*
 * Auxilliary functions
 */

void initDrugModule(int _withinHostTimestep, int _simulatorTimestep){
  DrugRegistry* registry;
  manager = ProteomeManager::getManager();
  registry = DrugRegistry::getRegistry();
  withinHostTimestep = _withinHostTimestep;
  simulatorTimestep  = _simulatorTimestep;
  Mutation* crt76 = manager->getMutation(string("CRT"), 76, 'T');
  Drug* s;
  // = new Drug("Sulfadoxine", "S", 0.1, 10*24*60); //Invented values
  //registry->addDrug(s);
  s = new Drug("Chloroquine", "CQ", 0.02, 45*24*60); //Based on Hoshen
  vector<Mutation*> crt76L = vector<Mutation*>();
  crt76L.push_back(crt76);
  s->addPDRule(crt76L, 204.0);
  s->addPDRule(vector<Mutation*>(), 68.0);
  s->parseProteomeInstances(manager);
  registry->addDrug(s);

}

/*
 * Drug and Dose code
 */

Dose::Dose(int _time, double _quantity) {
  time = _time;
  quantity = _quantity;
}

Dose::Dose(const Dose& _original) {
  time = _original.time;
  quantity = _original.quantity;
}

Dose::~Dose() {
}

Dose &Dose::operator=(const Dose &rhs) {
  this->time = rhs.time;
  this->quantity = rhs.quantity;
  return *this;
}

int Dose::operator==(const Dose &rhs) const {
  return (time==rhs.time) && (quantity==rhs.quantity);
}

int Dose::operator<(const Dose &rhs) const {
  return time<rhs.time;
}

ostream& operator<<(ostream& out, const Dose& dose) {
  out << dose.time << endl;
  out << dose.quantity << endl;
  return out;
}

Drug::Drug(const Drug &_original) {
  name = _original.name;
  abbreviation = _original.abbreviation;
  absorptionFactor = _original.absorptionFactor;
  halfLife = _original.halfLife;
  human = 0;
  doses = vector<Dose*>();
  _concentration = 0;
  requiredMutations = vector<vector<Mutation*> >();
  pdParameters = vector<double>();
  proteomePDParameters = map<int,double>();
}

Drug::Drug(string _name, string _abbreviation,
    double _absorptionFactor, double _halfLife) {
  name = _name;
  abbreviation = _abbreviation;
  absorptionFactor = _absorptionFactor;
  halfLife = _halfLife;
  human = 0;
  doses = vector<Dose*>();
  _concentration = 0;
  requiredMutations = vector<vector<Mutation*> >();
  pdParameters = vector<double>();
  proteomePDParameters = map<int,double>();
}

Drug::~Drug() {
  //We do not free human. That is OK.
  //We don't need to clear mutations as they are pointers to canonics, right?
  //delete requiredMutations;
  //delete pdParams;
  //delete proteomePDParams;
  vector<Dose*>::iterator it;
  for(it=doses.begin(); it!=doses.end(); it++) {
    delete (*it);
  }
}

istream& operator>>(istream& in, Drug& drug) {
  //Human is not saved. This is correct.
  ProteomeManager* manager = ProteomeManager::getManager();
  int numRules;
  int numCachedRules;
  int numDoses;
  int i;
  vector<Dose*>::const_iterator itD;
  vector<double>::const_iterator itPP;
  map<int,double>::const_iterator itPPD;
  in >> drug.name;
  in >> drug.abbreviation;
  in >> drug.absorptionFactor;
  in >> drug.halfLife; 
  in >> numRules;
  drug.requiredMutations = vector< vector<Mutation*> >();
  for (i=0; i<numRules; i++) {
    vector <Mutation*> rule = vector<Mutation*>();
    int numMuts;
    int j;
    in >> numMuts;
    for (j=0; j<numMuts; j++) {
      string proteinName;
      int position;
      char allele;
      in >> proteinName;
      in >> position;
      in >> allele;
      rule.push_back(
          manager->getMutation(proteinName, position, allele));
    }
    drug.requiredMutations.push_back(rule);
  }
  in >> numRules;
  drug.pdParameters = vector<double>();
  for (i=0; i<numRules; i++) {
    double parameter;
    in >> parameter;
    drug.pdParameters.push_back(parameter);
  }
  in >> numCachedRules;
  drug.proteomePDParameters = map<int, double>();
  for (i=0; i<numCachedRules; i++) {
    int proteomeID;
    double parameter;
    in >> proteomeID;
    in >> parameter;
    drug.proteomePDParameters[proteomeID] = parameter;
  }
  in >> drug._concentration;
  in >> drug._nextConcentration;
  in >> numDoses;
  drug.doses = vector<Dose*>();
  for (i=0; i<numDoses; i++) {
    int time;
    double quantity;
    in >> time;
    in >> quantity;
    drug.doses.push_back(new Dose(time, quantity));
  }
  return in;
}

ostream& operator<<(ostream& out, const Drug& drug) {
  //Human is not saved. This is correct.
  vector<Dose*>::const_iterator itD;
  vector< vector<Mutation*> >::const_iterator itRM;
  vector<double>::const_iterator itPP;
  map<int,double>::const_iterator itPPD;
  out << drug.name << endl;
  out << drug.abbreviation << endl;
  out << drug.absorptionFactor << endl;
  out << drug.halfLife << endl; 
  out << drug.requiredMutations.size() << endl;
  for (itRM=drug.requiredMutations.begin(); itRM!=drug.requiredMutations.end(); itRM++) {
    vector<Mutation*>::const_iterator itM;
    out << (*itRM).size() << endl;
    for (itM=(*itRM).begin(); itM!=(*itRM).end(); itM++) {
      out << (*itM)->getProteinName() << endl;
      out << (*itM)->getPosition() << endl;
      out << (*itM)->getAllele() << endl;
    }
  }
  out << drug.pdParameters.size() << endl;
  for (itPP=drug.pdParameters.begin(); itPP!=drug.pdParameters.end(); itPP++) {
    out << *itPP << endl ;
  }
  out << drug.proteomePDParameters.size() << endl;
  for (itPPD=drug.proteomePDParameters.begin(); itPPD!=drug.proteomePDParameters.end(); itPPD++) {
    out << (*itPPD).first << endl ;
    out << (*itPPD).second << endl ;
  }
  out << drug._concentration << endl;
  out << drug._nextConcentration << endl;
  out << drug.doses.size() << endl;
  for (itD=drug.doses.begin(); itD!=drug.doses.end(); itD++) {
    out << **itD;
  }
  return out;
}

Drug& Drug::operator=(const Drug &rhs) {
  this->name = rhs.name;
  this->abbreviation = rhs.abbreviation;
  this->requiredMutations = rhs.requiredMutations;
  this->pdParameters = rhs.pdParameters;
  this->proteomePDParameters = rhs.proteomePDParameters;

  this->human = human; //Maybe human should be copied? Probably not
  this->doses = doses; //Maybe list should be (shallow) copied? not know...
  this->_concentration = _concentration;
  this->_nextConcentration = _nextConcentration;
  return *this;
}

int Drug::operator==(const Drug &rhs) const {
  return abbreviation==rhs.abbreviation;
}

int Drug::operator<(const Drug &rhs) const {
  return abbreviation<rhs.abbreviation;
}

void Drug::setConcentration(double concentration) {
  _concentration = concentration;
  _nextConcentration = calculateDecay(withinHostTimestep);
}

void Drug::addConcentration(double concentration) {
  _concentration += concentration;
  _nextConcentration = calculateDecay(withinHostTimestep);
}

double Drug::calculateDrugFactor(Infection* infection) const {
  //Returning an average of 2 points
  //The wackiness below is because of const requirements (which operator[] on maps does not supply)
  double param = (*(proteomePDParameters.find(infection->getProteome()->getProteomeID()))).second;
  double startFactor = 3.8/(1+param/_concentration);
  double endFactor = 3.8/(1+param/_nextConcentration);
  return (startFactor + endFactor)/2;
}

double Drug::calculateDecay(int time) const {
  //k = log(2)/halfLife
  return _concentration * exp(-time*log(2.0)/halfLife);
}

void Drug::decay() {
  _concentration = _nextConcentration;
  _nextConcentration = calculateDecay(withinHostTimestep);
}


Drug Drug::use(Human* human) {
  // To be removed?
  Drug usedDrug = Drug(*this);
  human = human;
  return usedDrug;
}

void Drug::addPDRule(vector<Mutation*> ruleRequiredMutations, double pdFactor) {
  requiredMutations.push_back(ruleRequiredMutations);
  pdParameters.push_back(pdFactor);
}

void Drug::parseProteomeInstances(ProteomeManager* manager) {
  vector<ProteomeInstance*> instances = manager->getInstances();
  vector<ProteomeInstance*>::const_iterator it;
  int numRules = requiredMutations.size();
  for (it=instances.begin(); it !=instances.end(); it++) {
    //cerr << " Here goes instance";
    for(int rule=0; rule<numRules; rule++) {
      if ((*it)->hasMutations(requiredMutations[rule])) {
        proteomePDParameters[(*it)->getProteomeID()] = pdParameters[rule];
        //cerr << " rule: " << rule << "\n";
        break;
      }
    }
  }

}

/*
 * DrugProxy code
 */

DrugProxy::DrugProxy(Human* _human) {
  human = _human;
  registry = DrugRegistry::getRegistry();
}

void DrugProxy::medicate(string _drugAbbrev, double _qty, int _time) throw(int) {
  /* We ignore time for now (as it is only relevant for ACTs).
   *   As such, no doses are created, but concentration is updated.
   */
  //cerr << "Medicating with: " << _drugAbbrev << " " << _qty << "\n";
  list<Drug*>* drugs = human->getDrugs();
  Drug* myDrug = 0;
  list<Drug*>::iterator it;
  for (it=drugs->begin(); it!=drugs->end(); it++) {
    if ((*it)->getAbbreviation() == _drugAbbrev) {
      myDrug = (*it);
      break;
    }
  }
  if (myDrug==0) {
    myDrug = registry->getDrug(_drugAbbrev);
    drugs->push_back(myDrug);
  }
  myDrug->addConcentration(_qty*myDrug->getAbsorptionFactor()/human->getWeight());
}

double DrugProxy::calculateDrugsFactor(Infection* _infection) {
  //We will choose for now the smallest (ie, most impact)
  list<Drug*>::const_iterator it;
  double factor = 1; //no effect
  for (it=human->getDrugs()->begin(); it!=human->getDrugs()->end(); it++) {
    double drugFactor;
    drugFactor = (*it)->calculateDrugFactor(_infection);
    if (drugFactor< factor) {
      factor = drugFactor;
    }
  }
  return factor;
}

void DrugProxy::decayDrugs() {
  list<Drug*>::iterator it;
  for (it=human->getDrugs()->begin(); it!=human->getDrugs()->end(); it++) {
    (*it)->decay();
  }
}

ostream& operator<<(ostream& out, const DrugProxy& proxy) {
  //This will be needed in the future
  return out;
}

istream& operator>>(istream& in, DrugProxy& proxy) {
  //This will be needed in the future
  return in;
}


/*
 * DrugRegistry code
 */

DrugRegistry* DrugRegistry::instance = 0;
vector<Drug*> DrugRegistry::drugs = vector<Drug*>(); 

DrugRegistry::DrugRegistry() {
}

ostream& operator<<(ostream& out, const DrugRegistry& registry) {
  vector<Drug*>::const_iterator it;
  out << string("Drugs available:\n");
  for (it=registry.drugs.begin(); it!=registry.drugs.end(); it++) {
    out << string("  ") << **it;
    out << string("\n") ;
  }
  return out;
}

DrugRegistry* DrugRegistry::getRegistry() {
  if (instance == 0) {
    instance = new DrugRegistry;
  }
  return instance;
}

void DrugRegistry::addDrug(Drug* _drug) throw(int) {
  // In this case the good behaviour is HAVING an exception
  //   (We don't want to add a drug if it already exists)
  try {
    delete getDrug(_drug->getAbbreviation());
  }
  catch (int i) {
    drugs.push_back(_drug);
    return;
  }
  throw(1); //Element exists, we throw
}

Drug* DrugRegistry::getDrug(string _abbreviation) throw(int) {
  vector<Drug*>::iterator i;
  for(i=drugs.begin(); i!=drugs.end(); ++i) {
    if ((**i).getAbbreviation() == _abbreviation) {
        Drug* myDrug = *i;
        Drug* cloneDrug  = new Drug(*myDrug);
        return cloneDrug;
    }
  }
  throw(1);
}


