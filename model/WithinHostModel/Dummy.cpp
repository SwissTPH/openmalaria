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

// -----  Initialization  -----

DummyWithinHostModel::DummyWithinHostModel() :
    WithinHostModel(), _cumulativeh(0.0), _cumulativeY(0.0), _cumulativeYlag(0.0), _SPattenuationt(0), _MOI(0), patentInfections(0)
{
  W_GAUSS(0, sigma_i);	// FIXME: random call to keep these in sync
}

DummyWithinHostModel::DummyWithinHostModel(istream& in) :
    WithinHostModel(in)
{
  in >> _MOI; 
  in >> patentInfections; 
  in >> _SPattenuationt;
  in >> cumulativeY;
  in >> cumulativeh;
  in >> timeStepMaxDensity;
  in >> _cumulativeh;
  in >> _cumulativeY;
  in >> _cumulativeYlag;

  if ( _MOI <  0) {
    cerr << "Error reading checkpoint" << endl;
    exit(-3);
  }

  for(int i=0;i<_MOI;++i) {
    infections.push_back(DummyInfection(in));
  }

  if (Global::modelVersion & INCLUDES_PK_PD) {
    _proxy.read (in);
  }
}

DummyWithinHostModel::~DummyWithinHostModel() {
  clearAllInfections();
  if (Global::modelVersion & INCLUDES_PK_PD) {
    _proxy.destroy();
  }
}

// -----  Update function, called each step  -----

void DummyWithinHostModel::update (double age) {
  _proxy.setWeight (120.0 * wtprop[TransmissionModel::getAgeGroup(age)]);
  if (Global::modelVersion & INCLUDES_PK_PD) {
    std::list<DummyInfection>::iterator i;
    for(i=infections.begin(); i != infections.end(); i++){
      i->multiplyDensity(exp(-_proxy.calculateDrugsFactor(i->getProteome())));
    }
    _proxy.decayDrugs();
  }
}


// -----  Simple infection adders/removers  -----

void DummyWithinHostModel::newInfection(){
  if (_MOI <= 20) {
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
  _proxy.medicate(drugName, qty, time);
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
  human.setTotalDensity(0.0);
  human.setTimeStepMaxDensity(0.0);
  if (_cumulativeInfections >  0) {
    cumulativeh=_cumulativeh;
    cumulativeY=_cumulativeY;
    std::list<DummyInfection>::iterator i;
    for(i=infections.begin(); i!=infections.end(); i++){
      //std::cout<<"uis: "<<infData->duration<<std::endl;
      timeStepMaxDensity=human.getTimeStepMaxDensity();
      
      i->determineWithinHostDensity();
      timeStepMaxDensity=std::max((double)i->getDensity(), timeStepMaxDensity);
      human.setTimeStepMaxDensity(timeStepMaxDensity);
        
      human.setTotalDensity(human.getTotalDensity()+i->getDensity());
      //Compute the proportion of parasites remaining after innate blood stage effect
      if (i->getDensity() > Human::detectionlimit) {
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
}


// -----  Data checkpointing  -----

void DummyWithinHostModel::write(ostream& out) const {
  out << _cumulativeInfections << endl; 
  out << _pTransToMosq << endl;  
  out << _MOI << endl; 
  out << patentInfections << endl; 
  out << _SPattenuationt << endl;
  out << cumulativeY << endl;
  out << cumulativeh << endl;
  out << timeStepMaxDensity << endl;
  out << _cumulativeh << endl;
  out << _cumulativeY << endl;
  out << _cumulativeYlag << endl;

  for(std::list<DummyInfection>::const_iterator iter=infections.begin(); iter != infections.end(); iter++)
    iter->write (out);
  
  if (Global::modelVersion & INCLUDES_PK_PD) {
    _proxy.write (out);
  }
}
