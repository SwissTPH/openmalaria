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
#include "WithinHostModel/Empirical.h"
#include "simulation.h"
#include "TransmissionModel.h"	// getAgeGroup() - is in a funny place
#include "summary.h"
#include "inputData.h"

using namespace std;

const int EmpiricalWithinHostModel::MAX_INFECTIONS = 21;

// -----  Initialization  -----

EmpiricalWithinHostModel::EmpiricalWithinHostModel() :
    WithinHostModel(), drugProxy(DrugModel::createDrugModel ()),
    _MOI(0), patentInfections(0)
{
}

EmpiricalWithinHostModel::EmpiricalWithinHostModel(istream& in) :
    WithinHostModel(in), drugProxy(DrugModel::createDrugModel (in))
{
  in >> _MOI; 
  in >> patentInfections; 
  
  if (_MOI < 0 || _MOI > MAX_INFECTIONS)
    throw checkpoint_error ("_MOI");
  
  for(int i=0;i<_MOI;++i)
    infections.push_back(EmpiricalInfection(in));
}

EmpiricalWithinHostModel::~EmpiricalWithinHostModel() {
  clearAllInfections();
  delete drugProxy;
}

// -----  Update function, called each step  -----

void EmpiricalWithinHostModel::update () {
  std::list<EmpiricalInfection>::iterator i;
  for(i=infections.begin(); i != infections.end(); i++){
    i->multiplyDensity(exp(-drugProxy->getDrugFactor(i->getProteome())));
  }
  drugProxy->decayDrugs();
}


// -----  Simple infection adders/removers  -----

void EmpiricalWithinHostModel::newInfection(){
  if (_MOI < MAX_INFECTIONS) {
    _cumulativeInfections++;
    infections.push_back(EmpiricalInfection(Simulation::simulationTime, 1));
    _MOI++;
  }
}

void EmpiricalWithinHostModel::clearOldInfections(){
  std::list<EmpiricalInfection>::iterator iter=infections.begin();
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

void EmpiricalWithinHostModel::clearAllInfections(){
  std::list<EmpiricalInfection>::iterator i;
  for(i=infections.begin(); i != infections.end(); i++){
    i->destroy();
  }
  infections.clear();
  _MOI=0;
}


// -----  medicate drugs -----

void EmpiricalWithinHostModel::medicate(string drugName, double qty, int time, double age) {
  drugProxy->medicate(drugName, qty, time, 120.0 * wtprop[TransmissionModel::getAgeGroup(age)]);
}


// -----  Density calculations  -----

void EmpiricalWithinHostModel::calculateDensities(Human& human) {
  
  patentInfections = 0;
  totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  std::list<EmpiricalInfection>::iterator i;
  for(i=infections.begin(); i!=infections.end(); i++){
    i->determineWithinHostDensity();
    timeStepMaxDensity=std::max((double)i->getDensity(), timeStepMaxDensity);
    
    totalDensity += i->getDensity();
    //Compute the proportion of parasites remaining after innate blood stage effect
    if (i->getDensity() > detectionLimit) {
      patentInfections++;
    }
    _pTransToMosq = human.infectiousness();
  }
}

// -----  Summarize  -----

void EmpiricalWithinHostModel::summarize(double age) {
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

void EmpiricalWithinHostModel::write(ostream& out) const {
  out << _cumulativeInfections << endl; 
  out << _pTransToMosq << endl;  
  out << totalDensity << endl;
  out << timeStepMaxDensity << endl;
  
  drugProxy->write (out);
  
  out << _MOI << endl; 
  out << patentInfections << endl; 

  for(std::list<EmpiricalInfection>::const_iterator iter=infections.begin(); iter != infections.end(); iter++)
    iter->write (out);
}

//Immunity?

void EmpiricalWithinHostModel::updateImmuneStatus() {
}

void EmpiricalWithinHostModel::immunityPenalisation() {
}
