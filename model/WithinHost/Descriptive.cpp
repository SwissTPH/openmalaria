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
#include "WithinHost/Descriptive.h"
#include "Simulation.h"
#include "intervention.h"
#include "summary.h"

using namespace std;

const int DescriptiveWithinHostModel::MAX_INFECTIONS = 21;


// -----  Initialization  -----

DescriptiveWithinHostModel::DescriptiveWithinHostModel() :
    WithinHostModel(), _MOI(0),
    patentInfections(0)
{
  _innateImmunity = gsl::rngGauss(0, sigma_i);
}

DescriptiveWithinHostModel::DescriptiveWithinHostModel(istream& in) :
    WithinHostModel(in)
{
  readDescriptiveWHM (in);
  
  for(int i=0;i<_MOI;++i)
    infections.push_back(new DescriptiveInfection(in));
}

DescriptiveWithinHostModel::DescriptiveWithinHostModel(istream& in, bool) :
    WithinHostModel(in)
{
  readDescriptiveWHM (in);
}

DescriptiveWithinHostModel::~DescriptiveWithinHostModel() {
  clearAllInfections();
}


// -----  Data checkpointing  -----

void DescriptiveWithinHostModel::write(ostream& out) const {
  WithinHostModel::write (out);
  out << _MOI << endl;
  out << patentInfections << endl;
  out << _innateImmunity << endl;
  
  for(std::list<DescriptiveInfection*>::const_iterator iter=infections.begin(); iter != infections.end(); iter++)
    (*iter)->write (out);
}

void DescriptiveWithinHostModel::readDescriptiveWHM (istream& in) {
  in >> _MOI;
  in >> patentInfections; 
  in >> _innateImmunity;
  
  if (_MOI < 0 || _MOI > MAX_INFECTIONS)
    throw checkpoint_error ("_MOI");
}


// -----  Update function, called each step  -----

void DescriptiveWithinHostModel::update () {}


// -----  Simple infection adders/removers  -----

void DescriptiveWithinHostModel::newInfection(){
  if (_MOI < MAX_INFECTIONS) {
    infections.push_back(new DescriptiveInfection(Simulation::simulationTime));
    _MOI++;
  }
}

void DescriptiveWithinHostModel::clearAllInfections(){
  std::list<DescriptiveInfection*>::iterator i;
  for(i=infections.begin(); i != infections.end(); i++){
    delete *i;
  }
  infections.clear();
  _MOI=0;
}


// -----  Density calculations  -----

void DescriptiveWithinHostModel::calculateDensities(double ageInYears, double BSVEfficacy) {
  updateImmuneStatus ();	// inout(_cumulativeh,_cumulativeY)
  std::list<DescriptiveInfection*>::iterator iter=infections.begin();
  while(iter != infections.end()){
    if ((*iter)->expired()) {
      delete *iter;
      iter=infections.erase(iter);
      _MOI--;
    }
    else{
      iter++;
    }
  }//TODO cleanup
  
  patentInfections = 0;
  totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  
  // Values of _cumulativeh/Y at beginning of step
  // (values are adjusted for each infection)
  double cumulativeh=_cumulativeh;
  double cumulativeY=_cumulativeY;
  
  // IPTi SP dose clears infections at the time that blood-stage parasites appear     
  SPAction();
  
  for(iter=infections.begin(); iter!=infections.end(); iter++){
    //std::cout<<"uis: "<<infData->duration<<std::endl;
    // MAX_DENS_BUG: should be: infStepMaxDens = 0.0;
    double infStepMaxDens = timeStepMaxDensity;
    (*iter)->determineDensities(Simulation::simulationTime, ageInYears, cumulativeh, cumulativeY, infStepMaxDens, exp(-_innateImmunity), BSVEfficacy);
    
    IPTattenuateAsexualDensity (*iter);
    
    // MAX_DENS_BUG: should be: timeStepMaxDensity = std::max(infStepMaxDens, timeStepMaxDensity);
    timeStepMaxDensity = infStepMaxDens;
    
    totalDensity += (*iter)->getDensity();
    //Compute the proportion of parasites remaining after innate blood stage effect
    if ((*iter)->getDensity() > detectionLimit) {
      patentInfections++;
    }
    if ((*iter)->getStartDate() == (Simulation::simulationTime-1)) {
      _cumulativeh++;
    }
    (*iter)->determineDensityFinal ();
    _cumulativeY += Global::interval*(*iter)->getDensity();
  }
  
  IPTattenuateAsexualMinTotalDensity();
}

// -----  Summarize  -----

// TODO: can summarize move to WithinHostModel ?
void DescriptiveWithinHostModel::summarize(double age) {
  if (_MOI > 0) {
    Simulation::gMainSummary->addToInfectedHost(age,1);
    Simulation::gMainSummary->addToTotalInfections(age, _MOI);
    // TODO: patentInfections doesn't really need to be calculated each timestep and stored! calculate it here instead.
    Simulation::gMainSummary->addToTotalPatentInfections(age, patentInfections);
  }
  if (parasiteDensityDetectible()) {
    Simulation::gMainSummary->addToPatentHost(age, 1);
    Simulation::gMainSummary->addToSumLogDensity(age, log(totalDensity));
  }
}
