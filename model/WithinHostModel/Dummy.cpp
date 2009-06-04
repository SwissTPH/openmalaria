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

#include "GSLWrapper.h"
#include "human.h"
#include "WithinHostModel/Dummy.h"
#include "simulation.h"
#include "TransmissionModel.h"	// getAgeGroup() - is in a funny place
#include "summary.h"
#include "inputData.h"

using namespace std;

const int DummyWithinHostModel::MAX_INFECTIONS = 20;

// -----  Initialization  -----

DummyWithinHostModel::DummyWithinHostModel() :
    WithinHostModel(), drugProxy(DrugModel::createDrugModel ()),
    _cumulativeh(0.0), _cumulativeY(0.0), _cumulativeYlag(0.0),
    _MOI(0), patentInfections(0)
{
  W_GAUSS(0, sigma_i);	// FIXME: random call to keep these in sync
}

DummyWithinHostModel::DummyWithinHostModel(istream& in) :
    WithinHostModel(in), drugProxy(DrugModel::createDrugModel (in))
{
  in >> _MOI; 
  in >> patentInfections; 
  in >> cumulativeY;
  in >> cumulativeh;
  in >> _cumulativeh;
  in >> _cumulativeY;
  in >> _cumulativeYlag;
  
  if (_MOI < 0 || _MOI > MAX_INFECTIONS)
    throw checkpoint_error ("_MOI");
  
  for(int i=0;i<_MOI;++i)
    infections.push_back(DummyInfection(in));
}

DummyWithinHostModel::~DummyWithinHostModel() {
  clearAllInfections();
  delete drugProxy;
}

// -----  Update function, called each step  -----

void DummyWithinHostModel::update (double age) {
  drugProxy->setWeight (120.0 * wtprop[TransmissionModel::getAgeGroup(age)]);
  std::list<DummyInfection>::iterator i;
  for(i=infections.begin(); i != infections.end(); i++){
    i->multiplyDensity(exp(-drugProxy->getDrugFactor(i->getProteome())));
  }
  drugProxy->decayDrugs();
}


// -----  Simple infection adders/removers  -----

void DummyWithinHostModel::newInfection(){
  if (_MOI <= MAX_INFECTIONS) {
    _cumulativeInfections++;
    infections.push_back(DummyInfection(Simulation::simulationTime));
    _MOI++;
  }
}

void DummyWithinHostModel::clearOldInfections(){
  std::list<DummyInfection>::iterator iter=infections.begin();
  while(iter != infections.end()){
    int enddate=iter->getEndDate();
    if (Simulation::simulationTime >= enddate) {
      iter->destroy();
      iter=infections.erase(iter);
      _MOI--;
    }
    else{
      iter++;
    }
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

void DummyWithinHostModel::medicate(string drugName, double qty, int time) {
  drugProxy->medicate(drugName, qty, time);
}


// -----  immunity  -----

void DummyWithinHostModel::updateImmuneStatus(){
  if (immEffectorRemain < 1){
    _cumulativeh*=immEffectorRemain;
    _cumulativeY*=immEffectorRemain;
  }
  if (asexImmRemain < 1){
    _cumulativeh*=asexImmRemain/
        (1+(_cumulativeh*(1-asexImmRemain)/Infection::cumulativeHstar));
    _cumulativeY*=asexImmRemain/
        (1+(_cumulativeY*(1-asexImmRemain)/Infection::cumulativeYstar));
  }
}

void DummyWithinHostModel::immunityPenalisation() {
  _cumulativeY=(double)_cumulativeYlag-(immPenalty_22*(_cumulativeY-_cumulativeYlag));
  if (_cumulativeY <  0) {
    _cumulativeY=0.0;
  }
}


// -----  Density calculations  -----

void DummyWithinHostModel::calculateDensities(Human& human) {
  _cumulativeYlag = _cumulativeY;
  
  _pTransToMosq = 0.0;
  patentInfections = 0;
  totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  if (_cumulativeInfections >  0) {
    cumulativeh=_cumulativeh;
    cumulativeY=_cumulativeY;
    std::list<DummyInfection>::iterator i;
    for(i=infections.begin(); i!=infections.end(); i++){
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
  }
  _pTransToMosq = human.infectiousness();
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


// -----  Data checkpointing  -----

void DummyWithinHostModel::write(ostream& out) const {
  out << _cumulativeInfections << endl; 
  out << _pTransToMosq << endl;  
  out << totalDensity << endl;
  out << timeStepMaxDensity << endl;
  
  drugProxy->write (out);
  
  out << _MOI << endl; 
  out << patentInfections << endl; 
  out << cumulativeY << endl;
  out << cumulativeh << endl;
  out << _cumulativeh << endl;
  out << _cumulativeY << endl;
  out << _cumulativeYlag << endl;

  for(std::list<DummyInfection>::const_iterator iter=infections.begin(); iter != infections.end(); iter++)
    iter->write (out);
}
