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

#include "WithinHost/Dummy.h"
#include "inputData.h"
#include "util/errors.hpp"

using namespace std;

namespace OM { namespace WithinHost {
    // -----  Initialization  -----

DummyWithinHostModel::DummyWithinHostModel() :
    WithinHostModel(), pkpdModel(PkPd::PkPdModel::createPkPdModel ()),
    patentInfections(0)
{}

DummyWithinHostModel::~DummyWithinHostModel() {
  delete pkpdModel;
}


// -----  Simple infection adders/removers  -----

void DummyWithinHostModel::newInfection(){
  if (_MOI < MAX_INFECTIONS) {
        infections.push_back(DummyInfection());
    _MOI++;
  }
}

void DummyWithinHostModel::clearAllInfections(){
  infections.clear();
  _MOI=0;
}


// -----  medicate drugs -----

void DummyWithinHostModel::medicate(string drugName, double qty, int time, double age) {
  pkpdModel->medicate(drugName, qty, time, age);
}


// -----  Density calculations  -----

void DummyWithinHostModel::calculateDensities(double ageInYears, double BSVEfficacy) {
  updateImmuneStatus ();	// inout(_cumulativeh,_cumulativeY)
  
  patentInfections = 0;
  totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  
  for(std::list<DummyInfection>::iterator i=infections.begin(); i!=infections.end();) {
    if (Global::simulationTime >= i->getEndDate()) {
      i=infections.erase(i);
      _MOI--;
      continue;
    }
    
    double survivalFactor = (1.0-BSVEfficacy) * _innateImmSurvFact;
    survivalFactor *= pkpdModel->getDrugFactor(i->getProteome(), ageInYears);
    survivalFactor *= i->immunitySurvivalFactor(ageInYears, _cumulativeh, _cumulativeY);
    i->multiplyDensity(survivalFactor);
    i->determineWithinHostDensity();
    timeStepMaxDensity=std::max(i->getDensity(), timeStepMaxDensity);
    
    totalDensity += i->getDensity();
    //Compute the proportion of parasites remaining after innate blood stage effect
    if (i->getDensity() > detectionLimit) {
	patentInfections++;
    }
    if (i->getStartDate() == (Global::simulationTime-1)) {
	_cumulativeh++;
    }
    _cumulativeY += Global::interval*i->getDensity();
    
    ++i;
  }
  
  pkpdModel->decayDrugs();
}


// -----  Summarize  -----

int DummyWithinHostModel::countInfections (int& patentInfections) {
  if (infections.empty()) return 0;
  patentInfections = 0;
  for (std::list<DummyInfection>::iterator iter=infections.begin();
       iter != infections.end(); ++iter){
    if (iter->getDensity() > detectionLimit)
      patentInfections++;
  }
  return infections.size();
}


void DummyWithinHostModel::checkpoint (istream& stream) {
    WithinHostModel::checkpoint (stream);
    (*pkpdModel) & stream;
    patentInfections & stream;
    infections & stream;
    if (int(infections.size()) != _MOI)
	throw util::checkpoint_error ("_MOI mismatch");
}
void DummyWithinHostModel::checkpoint (ostream& stream) {
    WithinHostModel::checkpoint (stream);
    (*pkpdModel) & stream;
    patentInfections & stream;
    infections & stream;
}

} }