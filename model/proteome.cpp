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
#include <string>
#include <vector>

#include "proteome.h"

using namespace std;

/*
 * Support functions
 */

void initProteomeModule() {
  Protein* crt = new Protein(string("CRT"));
  ProteinPosition* pos = new ProteinPosition(crt, 76, 'K');
  Mutation* mutation = new Mutation(pos, 'T');
  ProteomeInstance* instance;
  instance = new ProteomeInstance();
  ProteomeManager::addInstance(instance);
  instance = new ProteomeInstance();
  instance->addMutation(mutation);
  ProteomeManager::addInstance(instance);
  //ProteomeManager* manager;
  //manager = ProteomeManager::getManager();
}


/*
 * Protein
 */

Protein::Protein(string _name) {
  name = _name;
  positions = vector<ProteinPosition*>();
  ProteomeManager::addProtein(this);
}

Protein::~Protein() {
  vector<ProteinPosition*>::const_iterator it;
  for(it=positions.begin(); it!=positions.end(); it++) {
    delete (*it);
  }
}

string Protein::getName() const {
  return name;
}

void Protein::addPosition(ProteinPosition* _position) {
  positions.push_back(_position);
}

Mutation* Protein::getMutation(int _position, char _allele) throw(int) {
  vector<ProteinPosition*>::const_iterator it;
  for (it=positions.begin(); it != positions.end(); it++) {
    if ((*it)->getPosition() == _position) {
      return (*it)->getMutation(_allele);
    }
  }
  throw(2); // Position not known
}

ostream& operator<<(ostream& out, const Protein& protein) {
  vector<ProteinPosition*>::const_iterator it;
  out << protein.name << endl;
  out << protein.positions.size() << endl;
  for(it=protein.positions.begin(); it != protein.positions.end(); it++) {
    vector<Mutation*>::const_iterator itM;
    ProteinPosition* position = (*it);
    out << position->getPosition() << endl;
    out << position->getWildType() << endl;
    out << position->getMutations()->size() << endl;
    for(itM=position->getMutations()->begin(); itM != position->getMutations()->end(); itM++) {
      out << (*itM)->getAllele() << endl;
    }
  }
  return out;
}

istream& operator>>(istream& in, Protein& protein) {
  int numPositions;
  int numMutations;
  int i;
  int j;
  in >> protein.name;
  in >> numPositions;
  protein.positions = vector<ProteinPosition*>();
  for(i=0; i<numPositions; i++) {
    int pos;
    char wildType;
    in >> pos;
    in >> wildType;
    //Object will add itself to protein.positions
    ProteinPosition* position = new ProteinPosition(&protein, pos, wildType);
    in >> numMutations;
    for(j=0; j<numMutations; j++) {
      char allele;
      in >> allele;
      new Mutation(position, allele); //it will hook itself to position
    }
  }
  return in;
}


/*
 * ProteinPosition
 */

ProteinPosition::ProteinPosition(Protein* _protein, int _position,
    char _wildType) {
  protein = _protein;
  protein->addPosition(this);
  position = _position;
  wildType = _wildType;
  mutations = vector<Mutation*>();
}

ProteinPosition::~ProteinPosition() {
  vector<Mutation*>::const_iterator it;
  for(it=mutations.begin(); it!=mutations.end(); it++) {
    delete (*it);
  }
}

Protein* ProteinPosition::getProtein() {
  return protein;
}

string ProteinPosition::getProteinName() {
  return protein->getName();
}

ostream& operator<<(ostream& out, const ProteinPosition& position) {
  return out;
}

istream& operator>>(istream& in, ProteinPosition& position) {
  return in;
}

void ProteinPosition::addMutation(Mutation* _mutation) {
  mutations.push_back(_mutation);
}

int ProteinPosition::getPosition() const {
  return position;
}

/*
 * Mutation
 */

Mutation* ProteinPosition::getMutation(char _allele) throw(int) {
  vector<Mutation*>::const_iterator it;
  for (it=mutations.begin(); it != mutations.end(); it++) {
    if ((*it)->getAllele() == _allele) {
      return *it;
    }
  }
  throw(3); // Allele not known
}

Mutation::Mutation(ProteinPosition* _position, char _allele) {
  allele = _allele;
  position = _position;
  position->addMutation(this);
}

Mutation::~Mutation() {
}

ostream& operator<<(ostream& out, const Mutation& mutation) {
  return out;
}

istream& operator>>(istream& in, Mutation& mutation) {
  return in;
}


char Mutation::getAllele() const {
  return allele;
}

int Mutation::getPosition() const {
  return position->getPosition();
}

string Mutation::getProteinName() const {
  return position->getProteinName();
}

int Mutation::operator==(const Mutation& rhs) const {
  if (allele==rhs.allele) {
    if (position->getPosition() == rhs.position->getPosition()) {
      if (position->getProtein()->getName() == rhs.position->getProtein()->getName()) {
        return 1;
      }
    }
  }
  return 0;
}

/*
 * ProteomeInstance
 */

int ProteomeInstance::currentID = 0;

ProteomeInstance::ProteomeInstance(){
  proteomeID = currentID++;
  mutations = vector<Mutation*>();
}

ProteomeInstance::~ProteomeInstance(){
  //the content doesnt need to be deleted, that is the reponsability of delete proteins
  mutations.clear();
}

ostream& operator<<(ostream& out, const ProteomeInstance& instance) {
  vector<Mutation*>::const_iterator it;
  out << instance.proteomeID << endl;
  out << instance.mutations.size() << endl;
  for(it=instance.mutations.begin(); it!=instance.mutations.end(); it++) {
    out << (*it)->getProteinName() << endl;
    out << (*it)->getPosition() << endl;
    out << (*it)->getAllele() << endl;
  }
  return out;
}

istream& operator>>(istream& in, ProteomeInstance& instance) {
  int numMuts;
  in >> instance.proteomeID;
  in >> numMuts;
  instance.mutations = vector<Mutation*>();
  for (int i=0; i<numMuts; i++) {
    string proteinName;
    int position;
    char allele;
    in >> proteinName;
    in >> position;
    in >> allele;
    instance.mutations.push_back(
      ProteomeManager::getManager()->getMutation(proteinName, position, allele));
  }
  if (instance.currentID<instance.proteomeID) {
    instance.currentID = instance.proteomeID;
  }
  return in;
}

void ProteomeInstance::addMutation(Mutation* _mutation) {
  mutations.push_back(_mutation);
}

int ProteomeInstance::getProteomeID() {
  return proteomeID;
}

bool ProteomeInstance::hasMutations(vector<Mutation*> _mutations) {
  vector<Mutation*>::const_iterator it;
  for(it=_mutations.begin(); it!=_mutations.end(); it++) {
    vector<Mutation*>::const_iterator it2;
    for(it2=mutations.begin(); it2!=mutations.end(); it2++) {
      if ((**it)==(**it2)) {
        goto findNext; //Yes, it is a goto, and a very good one.
      }
    }
    return false;
    findNext: continue;
  }
  return true;
}

/*
 * ProteomeManager
 */

ProteomeManager* ProteomeManager::manager = 0;
vector<ProteomeInstance*> ProteomeManager::instances = vector<ProteomeInstance*>();
vector<Protein*> ProteomeManager::proteins = vector<Protein*>();


ProteomeManager::ProteomeManager() {
  manager = this;
}

ProteomeManager::~ProteomeManager() {
}

ostream& operator<<(ostream& out, const ProteomeManager& manager) {
  vector<ProteomeInstance*>* instances = &ProteomeManager::instances;
  vector<Protein*>* proteins = &ProteomeManager::proteins;
  vector<ProteomeInstance*>::const_iterator iti;
  vector<Protein*>::const_iterator itp;

  out << proteins->size() << endl;
  for(itp=proteins->begin(); itp!=proteins->end(); itp++) {
    out << *(*itp);
  }
  out << instances->size() << endl;
  for(iti=instances->begin(); iti!=instances->end(); iti++) {
    out << *(*iti);
  }
  return out;
}

istream& operator>>(istream& in, ProteomeManager& manager) {
  //Proteins need to go before instances!
  int numInstances;
  int numProteins;
  int i;
  if (ProteomeManager::manager != 0) { //This is aways true...
    vector<ProteomeInstance*>::const_iterator iti ;
    vector<Protein*>::const_iterator itp;
    for(iti=ProteomeManager::instances.begin(); iti!=ProteomeManager::instances.end(); iti++) {
      delete (*iti);
    }
    for(itp=ProteomeManager::proteins.begin(); itp!=ProteomeManager::proteins.end(); itp++) {
      delete (*itp);
    }
    ProteomeManager::instances.clear();
    ProteomeManager::proteins.clear();
  }
  ProteomeManager::manager = &manager;
  ProteomeManager::instances = vector<ProteomeInstance*>();
  ProteomeManager::proteins = vector<Protein*>();
  in >> numProteins;
  for(i=0; i<numProteins; i++) {
    Protein* protein = new Protein("");
    in >> (*protein);
  }
  in >> numInstances;
  for(i=0; i<numInstances; i++) {
    ProteomeInstance* instance = new ProteomeInstance();
    in >> (*instance);
    ProteomeManager::instances.push_back(instance);
  }
  return in;
}


void ProteomeManager::addInstance(ProteomeInstance* _instance) {
  instances.push_back(_instance);
}

void ProteomeManager::addProtein(Protein* _protein) {
  proteins.push_back(_protein);
}

ProteomeManager* ProteomeManager::getManager() {
  if (manager==0) {
    ProteomeManager();
  }
  return manager;
}

Mutation* ProteomeManager::getMutation(string _proteinName, int _position, char _allele) throw(int) {
  vector<Protein*>::const_iterator itProt;
  for (itProt=proteins.begin(); itProt != proteins.end(); itProt++) {
    if ((*itProt)->getName() == _proteinName) {
      return (*itProt)->getMutation(_position, _allele);
    }
  }
  throw(1); // Name not known
}

vector<ProteomeInstance*> ProteomeManager::getInstances() const {
  return instances;
}

ProteomeInstance* ProteomeManager::getProteome(int proteome) const {
  return instances[proteome];
}

ProteomeInstance* ProteomeManager::getInfection() const {
  return *(instances.begin()); //To be changed
}

