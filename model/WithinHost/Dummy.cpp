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

#include "util/gsl.h"
#include "WithinHost/Dummy.h"
#include "Simulation.h"
#include "summary.h"
#include "inputData.h"

using namespace std;

// -----  Initialization  -----

DummyWithinHostModel::DummyWithinHostModel() :
    WithinHostModel(), drugProxy(DrugModel::createDrugModel ()),
    _cumulativeh(0.0), _cumulativeY(0.0), _cumulativeYlag(0.0),
    _MOI(0), patentInfections(0)
{}

DummyWithinHostModel::~DummyWithinHostModel() {
  clearAllInfections();
  delete drugProxy;
}

DummyWithinHostModel::DummyWithinHostModel(istream& in) :
    WithinHostModel(in), drugProxy(DrugModel::createDrugModel (in))
{
  in >> _MOI; 
  in >> patentInfections; 
  in >> _cumulativeh;
  in >> _cumulativeY;
  in >> _cumulativeYlag;
  
  if (_MOI < 0 || _MOI > MAX_INFECTIONS)
    throw checkpoint_error ("_MOI");
  
  for(int i=0;i<_MOI;++i)
    infections.push_back(DummyInfection(in));
}
void DummyWithinHostModel::write(ostream& out) const {
  WithinHostModel::write (out);
  drugProxy->write (out);
  
  out << _MOI << endl; 
  out << patentInfections << endl; 
  out << _cumulativeh << endl;
  out << _cumulativeY << endl;
  out << _cumulativeYlag << endl;
  
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
  drugProxy->medicate(drugName, qty, time, age, 120.0 * wtprop[getAgeGroup(age)]);
}


// -----  Density calculations  -----

void DummyWithinHostModel::calculateDensities(double ageInYears, double BSVEfficacy) {
  updateImmuneStatus ();	// inout(_cumulativeh,_cumulativeY)
  
  patentInfections = 0;
  totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  for(std::list<DummyInfection>::iterator i=infections.begin(); i!=infections.end(); i++) {
    if (Simulation::simulationTime >= i->getEndDate()) {
      i->destroy();
      i=infections.erase(i);
      _MOI--;
      continue;
    }
    
    i->multiplyDensity(drugProxy->getDrugFactor(i->getProteome()));
    i->determineWithinHostDensity();
    timeStepMaxDensity=std::max((double)i->getDensity(), timeStepMaxDensity);
    
    totalDensity += i->getDensity();
    //Compute the proportion of parasites remaining after innate blood stage effect
    if (i->getDensity() > detectionLimit) {
      patentInfections++;
    }
    if (i->getStartDate() == (Simulation::simulationTime-1)) {
      _cumulativeh++;
    }
    _cumulativeY += Global::interval*i->getDensity();
  }
  
  drugProxy->decayDrugs();
}

// -----  Summarize  -----

void DummyWithinHostModel::summarize(double age) {
  if (_MOI > 0) {
    Simulation::gMainSummary->addToInfectedHost(age,1);
    Simulation::gMainSummary->addToTotalInfections(age, _MOI);
    Simulation::gMainSummary->addToTotalPatentInfections(age, patentInfections);
  }
  if (parasiteDensityDetectible()) {
    Simulation::gMainSummary->addToPatentHost(age, 1);
    Simulation::gMainSummary->addToSumLogDensity(age, log(totalDensity));
  }
}
