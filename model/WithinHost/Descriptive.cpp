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

#include "WithinHost/Descriptive.h"
#include "intervention.h"

using namespace std;


// -----  Initialization  -----

DescriptiveWithinHostModel::DescriptiveWithinHostModel() :
    WithinHostModel()
{}

DescriptiveWithinHostModel::DescriptiveWithinHostModel(istream& in) :
    WithinHostModel(in)
{
  for(int i=0;i<_MOI;++i)
    infections.push_back(new DescriptiveInfection(in));
}

DescriptiveWithinHostModel::DescriptiveWithinHostModel(istream& in, bool) :
    WithinHostModel(in)
{}

DescriptiveWithinHostModel::~DescriptiveWithinHostModel() {
  clearAllInfections();
}


// -----  Data checkpointing  -----

void DescriptiveWithinHostModel::write(ostream& out) const {
  WithinHostModel::write (out);
  
  for(std::list<DescriptiveInfection*>::const_iterator iter=infections.begin(); iter != infections.end(); iter++)
    (*iter)->write (out);
}


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
  
  totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  
  // Values of _cumulativeh/Y at beginning of step
  // (values are adjusted for each infection)
  double cumulativeh=_cumulativeh;
  double cumulativeY=_cumulativeY;
  
  // IPTi SP dose clears infections at the time that blood-stage parasites appear     
  SPAction();
  
  for(iter=infections.begin(); iter!=infections.end(); iter++){
    // With option MAX_DENS_RESET this would be: infStepMaxDens = 0.0;
    // However, when using MAX_DENS_CORRECTION this is irrelevant.
    double infStepMaxDens = timeStepMaxDensity;
    (*iter)->determineDensities(ageInYears, cumulativeh, cumulativeY, infStepMaxDens, _innateImmSurvFact, BSVEfficacy);
    
    IPTattenuateAsexualDensity (*iter);
    
    if (Global::modelVersion & MAX_DENS_CORRECTION)
	infStepMaxDens = std::max(infStepMaxDens, timeStepMaxDensity);
    timeStepMaxDensity = infStepMaxDens;
    
    totalDensity += (*iter)->getDensity();
    if ((*iter)->getStartDate() == (Simulation::simulationTime-1)) {
      _cumulativeh++;
    }
    (*iter)->determineDensityFinal ();
    _cumulativeY += Global::interval*(*iter)->getDensity();
  }
  
  IPTattenuateAsexualMinTotalDensity();
}


// -----  Summarize  -----

int DescriptiveWithinHostModel::countInfections (int& patentInfections) {
  if (infections.empty()) return 0;
  patentInfections = 0;
  for (std::list<DescriptiveInfection*>::iterator iter=infections.begin();
       iter != infections.end(); ++iter){
    if ((*iter)->getDensity() > detectionLimit)
      patentInfections++;
  }
  return infections.size();
}
