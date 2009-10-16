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

#include <string>
#include <vector>

#include "proteome.h"
#include "util/gsl.h"

using namespace std;

/*
 * Support functions
 */

void initProteomeModule() {
  Protein* crt = new Protein(string("CRT"));
  ProteinPosition* pos = new ProteinPosition(crt, 76, 'K');
  Mutation* mutation = new Mutation(pos, 'T');
  ProteomeManager::addProtein(crt);
  
  ProteomeInstance::init (mutation);
}


/*
 * Protein
 */

Protein::Protein(string _name) {
  name = _name;
  positions = vector<ProteinPosition*>();
}

Protein::Protein(istream& in) {
  in >> name;
  positions = vector<ProteinPosition*>();
  size_t numPositions;
  in >> numPositions;
  Global::validateListSize (numPositions);
  for(size_t i=0; i<numPositions; i++) {
    //Object will add itself to this->positions
    new ProteinPosition (this, in);
  }
}

Protein::~Protein() {
  vector<ProteinPosition*>::const_iterator it;
  for(it=positions.begin(); it!=positions.end(); it++) {
    delete (*it);
  }
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

void Protein::write (ostream& out) {
  vector<ProteinPosition*>::const_iterator it;
  out << name << endl;
  out << positions.size() << endl;
  for(it=positions.begin(); it != positions.end(); it++) {
    (*it)->write(out);
  }
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
}
ProteinPosition::ProteinPosition(Protein* _protein, istream& in) {
  protein = _protein;
  protein->addPosition(this);
  in >> position;
  in >> wildType;
  size_t numMutations;
  in >> numMutations;
  Global::validateListSize (numMutations);
  for(size_t j=0; j<numMutations; j++) {
    char allele;
    in >> allele;
    new Mutation(this, allele); //it will hook itself to position
  }
}

ProteinPosition::~ProteinPosition() {
  vector<Mutation*>::const_iterator it;
  for(it=mutations.begin(); it!=mutations.end(); it++) {
    delete (*it);
  }
}

void ProteinPosition::write (ostream& out) {
  out << position << endl;
  out << wildType << endl;
  out << mutations.size() << endl;
  for(vector<Mutation*>::const_iterator itM=mutations.begin(); itM != mutations.end(); itM++) {
    out << (*itM)->getAllele() << endl;
  }
}

void ProteinPosition::addMutation(Mutation* _mutation) {
  mutations.push_back(_mutation);
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


// -----  ProteomeInstance static  -----

int ProteomeInstance::currentID;
vector<ProteomeInstance> ProteomeInstance::instances;

void ProteomeInstance::init (Mutation* mutation) {
  currentID = 0;
  // We store objects, not pointers - so this creates two objects, using the
  // default constructor. Then we don't need to free memory each object.
  // Note: can't add both at once, as the constructors must be called in the right order.
  instances.resize (1);
  instances.resize (2);
  instances[1].mutations.push_back(mutation);
}


const ProteomeInstance* ProteomeInstance::newInfection() {
  const double percRes = 0.0;		// proportion of resistant infections
  if (gsl::rngUniform() < percRes) {
    return &instances[1];
  } else {
    return &instances[0];
  }
}

/* Instances are created the same way on start-up, so no need to checkpoint
 * (for now at least).
 * Code to write:
  out << instances.size() << endl;
  for(vector<ProteomeInstance*>::const_iterator iti=instances.begin(); iti!=instances.end(); iti++) {
    (*iti)->write (out);
  }
 * Code to read:
  if (instances.size() > 0) {
    throw ...;
  }
  size_t numInstances;
  in >> numInstances;
  Global::validateListSize (numInstances);
  for(size_t i=0; i<numInstances; i++) {
    instances.push_back(new ProteomeInstance(in));
  }
 */


// -----  ProteomeInstance non-static  -----

ProteomeInstance::ProteomeInstance() {
  proteomeID = currentID++;
}

/*ProteomeInstance::ProteomeInstance (istream& in) {
  in >> proteomeID;
  size_t numMuts;
  in >> numMuts;
  Global::validateListSize (numMuts);
  mutations = vector<Mutation*>();
  for (size_t i=0; i<numMuts; i++) {
    string proteinName;
    int position;
    char allele;
    in >> proteinName;
    in >> position;
    in >> allele;
    mutations.push_back(
      ProteomeManager::getMutation(proteinName, position, allele));
  }
  if (currentID<proteomeID)
    currentID = proteomeID;
}*/

ProteomeInstance::~ProteomeInstance(){
  //the content doesnt need to be deleted, that is the reponsability of delete proteins
  mutations.clear();
}

/*void ProteomeInstance::write (ostream& out) {
  out << proteomeID << endl;
  out << mutations.size() << endl;
  for(vector<Mutation*>::const_iterator it=mutations.begin(); it!=mutations.end(); it++) {
    out << (*it)->getProteinName() << endl;
    out << (*it)->getPosition() << endl;
    out << (*it)->getAllele() << endl;
  }
}*/

bool ProteomeInstance::hasMutations(vector<Mutation*> _mutations) const {
  vector<Mutation*>::const_iterator it;
  for(it=_mutations.begin(); it!=_mutations.end(); it++) {
    list<Mutation*>::const_iterator it2;
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

vector<Protein*> ProteomeManager::proteins = vector<Protein*>();


void ProteomeManager::write (ostream& out) {
  out << proteins.size() << endl;
  for(vector<Protein*>::const_iterator itp=proteins.begin(); itp!=proteins.end(); itp++) {
    (*itp)->write (out);
  }
}

void ProteomeManager::read (istream& in) {
    for(vector<Protein*>::const_iterator itp=proteins.begin(); itp!=proteins.end(); itp++) {
      delete (*itp);
    }
    proteins.clear();
  
  size_t numProteins;
  in >> numProteins;
  Global::validateListSize (numProteins);
  for(size_t i=0; i<numProteins; i++) {
    proteins.push_back(new Protein(in));
  }
}


void ProteomeManager::addProtein(Protein* _protein) {
  proteins.push_back(_protein);
}

Mutation* ProteomeManager::getMutation(string _proteinName, int _position, char _allele) throw(int) {
  for (vector<Protein*>::const_iterator itProt = proteins.begin(); itProt != proteins.end(); itProt++) {
    if ((*itProt)->getName() == _proteinName) {
      return (*itProt)->getMutation(_position, _allele);
    }
  }
  throw(1); // Name not known
}

