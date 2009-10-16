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
#include "WithinHost/Empirical.h"
#include "Simulation.h"
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
EmpiricalWithinHostModel::~EmpiricalWithinHostModel() {
  clearAllInfections();
  delete drugProxy;
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
void EmpiricalWithinHostModel::write(ostream& out) const {
  writeWHM (out);
  drugProxy->write (out);
  
  out << _MOI << endl; 
  out << patentInfections << endl; 
  
  for(std::list<EmpiricalInfection>::const_iterator iter=infections.begin(); iter != infections.end(); iter++)
    iter->write (out);
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

void EmpiricalWithinHostModel::clearAllInfections(){
  infections.clear();
  _MOI=0;
}


// -----  medicate drugs -----

void EmpiricalWithinHostModel::medicate(string drugName, double qty, int time, double age) {
  drugProxy->medicate(drugName, qty, time, age, 120.0 * wtprop[getAgeGroup(age)]);
}


// -----  Density calculations  -----

void EmpiricalWithinHostModel::calculateDensities(double ageInYears, double BSVEfficacy) {
  patentInfections = 0;
  totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  std::list<EmpiricalInfection>::iterator i;
  for(i=infections.begin(); i!=infections.end();){
    if (i->updateDensity(Simulation::simulationTime)) {
      i = infections.erase(i);
      --_MOI;
      continue;
    }
    
    //FIXME: timeStepMaxDensity was meant to be the highest total density of any day within a 5-day interval wasn't it?
    // In which case, here it should be equal to totalDensity!
    timeStepMaxDensity=std::max((double)i->getDensity(), timeStepMaxDensity);
    totalDensity += i->getDensity();
    //Compute the proportion of parasites remaining after innate blood stage effect
    if (i->getDensity() > detectionLimit) {
      patentInfections++;
    }
    ++i;
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


//Immunity?

void EmpiricalWithinHostModel::updateImmuneStatus() {
}

void EmpiricalWithinHostModel::immunityPenalisation() {
}
