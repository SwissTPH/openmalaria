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

#include "PkPd/Proteome.h"
#include "util/gsl.h"

#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

namespace OM { namespace PkPd {
    
/*
 * Protein
 */

Protein::~Protein () {
    vector<ProteinPosition*>::const_iterator it;
    for (it=positions.begin(); it != positions.end(); it++)
	delete *it;
}

void Protein::addPosition(ProteinPosition* _position) {
  positions.push_back(_position);
}

Mutation* Protein::getMutation(int _position, char _allele) {
  vector<ProteinPosition*>::const_iterator it;
  for (it=positions.begin(); it != positions.end(); it++) {
    if ((*it)->getPosition() == _position) {
      return (*it)->getMutation(_allele);
    }
  }
  throw runtime_error ("Position not known");
}


/*
 * ProteinPosition
 */

ProteinPosition::ProteinPosition(Protein* _protein, int _position, char _wildType) {
  protein = _protein;
  protein->addPosition(this);
  position = _position;
  wildType = _wildType;
}

ProteinPosition::~ProteinPosition() {
  vector<Mutation*>::const_iterator it;
  for(it=mutations.begin(); it!=mutations.end(); it++) {
    delete (*it);
  }
}

void ProteinPosition::addMutation(Mutation* _mutation) {
  mutations.push_back(_mutation);
}


/*
 * Mutation
 */

Mutation* ProteinPosition::getMutation(char _allele) {
  vector<Mutation*>::const_iterator it;
  for (it=mutations.begin(); it != mutations.end(); it++) {
    if ((*it)->getAllele() == _allele) {
      return *it;
    }
  }
  throw runtime_error("Allele not known");
}

Mutation::Mutation(ProteinPosition* _position, char _allele) {
  allele = _allele;
  position = _position;
  position->addMutation(this);
}

Mutation::~Mutation() {
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
void ProteomeInstance::cleanup () {
    instances.resize (0);
}


const ProteomeInstance* ProteomeInstance::newInfection() {
  const double percRes = 0.0;		// proportion of resistant infections
  if (gsl::rngUniform() < percRes) {
    return &instances[1];
  } else {
    return &instances[0];
  }
}


// -----  ProteomeInstance non-static  -----

ProteomeInstance::ProteomeInstance() {
  proteomeID = currentID++;
}

ProteomeInstance::~ProteomeInstance(){
  //the content doesnt need to be deleted, that is the reponsability of delete proteins
  mutations.clear();
}

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


void ProteomeManager::init () {
    Protein* crt = new Protein(string("CRT"));
    ProteinPosition* pos = new ProteinPosition(crt, 76, 'K');
    Mutation* mutation = new Mutation(pos, 'T');
    ProteomeManager::addProtein(crt);
    
    ProteomeInstance::init (mutation);
}
void ProteomeManager::cleanup () {
    ProteomeInstance::cleanup ();
    
    for(vector<Protein*>::const_iterator itp=proteins.begin(); itp!=proteins.end(); itp++)
	delete *itp;
    proteins.resize (0);
}


void ProteomeManager::addProtein(Protein* _protein) {
  proteins.push_back(_protein);
}

Mutation* ProteomeManager::getMutation(string _proteinName, int _position, char _allele) {
  for (vector<Protein*>::const_iterator itProt = proteins.begin(); itProt != proteins.end(); itProt++) {
    if ((*itProt)->isNamed(_proteinName)) {
      return (*itProt)->getMutation(_position, _allele);
    }
  }
  throw runtime_error ("Name not known");
}

} }