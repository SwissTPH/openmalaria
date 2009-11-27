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
#include "Simulation.h"
#include "inputData.h"

using namespace std;

// -----  Initialization  -----

DummyWithinHostModel::DummyWithinHostModel() :
    WithinHostModel(), pkpdModel(PkPdModel::createPkPdModel ()),
    _cumulativeh(0.0), _cumulativeY(0.0), _cumulativeYlag(0.0),
    patentInfections(0)
{}

DummyWithinHostModel::~DummyWithinHostModel() {
  clearAllInfections();
  delete pkpdModel;
}

DummyWithinHostModel::DummyWithinHostModel(istream& in) :
    WithinHostModel(in), pkpdModel(PkPdModel::createPkPdModel (in))
{
  in >> patentInfections; 
  
  for(int i=0;i<_MOI;++i)
    infections.push_back(DummyInfection(in));
}
void DummyWithinHostModel::write(ostream& out) const {
  WithinHostModel::write (out);
  pkpdModel->write (out);
  
  out << patentInfections << endl; 
  
  for(std::list<DummyInfection>::const_iterator iter=infections.begin(); iter != infections.end(); iter++)
    iter->write (out);
}


// -----  Simple infection adders/removers  -----

void DummyWithinHostModel::newInfection(){
  if (_MOI < MAX_INFECTIONS) {
        infections.push_back(DummyInfection(Simulation::simulationTime));
    _MOI++;
  }
}

void DummyWithinHostModel::clearAllInfections(){
  std::list<DummyInfection>::iterator i;
  for(i=infections.begin(); i != infections.end(); i++){
    i->destroy();
  }
  infections.clear();
  _MOI=0;
}


// -----  medicate drugs -----

void DummyWithinHostModel::medicate(string drugName, double qty, int time, double age) {
  pkpdModel->medicate(drugName, qty, time, age, 120.0 * wtprop[getAgeGroup(age)]);
}


// -----  Density calculations  -----

void DummyWithinHostModel::calculateDensities(double ageInYears, double BSVEfficacy) {
  updateImmuneStatus ();	// inout(_cumulativeh,_cumulativeY)
  
  patentInfections = 0;
  totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  for(std::list<DummyInfection>::iterator i=infections.begin(); i!=infections.end(); ) {
    if (Simulation::simulationTime >= i->getEndDate()) {
      i=infections.erase(i);
      _MOI--;
      continue;
    }
    else {
	double survivalFactor = (1.0-BSVEfficacy) * _innateImmSurvFact;
	survivalFactor *= pkpdModel->getDrugFactor(i->getProteome());
	survivalFactor *= i->immunitySurvivalFactor(ageInYears, _cumulativeh, _cumulativeY);
	i->multiplyDensity(survivalFactor);
      i->determineWithinHostDensity();
      timeStepMaxDensity=std::max(i->getDensity(), timeStepMaxDensity);
    
      totalDensity += i->getDensity();
      //Compute the proportion of parasites remaining after innate blood stage effect
      if (i->getDensity() > detectionLimit) {
        patentInfections++;
      }
      if (i->getStartDate() == (Simulation::simulationTime-1)) {
        _cumulativeh++;
      }
      _cumulativeY += Global::interval*i->getDensity();
	  if (i == infections.end()) {
        break;
      }
      else {
        i++;
      }
    }
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
