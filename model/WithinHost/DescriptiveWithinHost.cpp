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
#include "util/ModelOptions.hpp"

using namespace std;

namespace OM { namespace WithinHost {
    
// -----  Initialization  -----

DescriptiveWithinHostModel::DescriptiveWithinHostModel() :
    WithinHostModel()
{}

DescriptiveWithinHostModel::~DescriptiveWithinHostModel() {
  clearAllInfections();
}


// -----  Simple infection adders/removers  -----

void DescriptiveWithinHostModel::newInfection(){
  if (_MOI < MAX_INFECTIONS) {
    infections.push_back(new DescriptiveInfection());
    _MOI++;
  }
}
void DescriptiveWithinHostModel::loadInfection(istream& stream){
    infections.push_back(new DescriptiveInfection(stream));
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
  updateImmuneStatus ();
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
    
    if (util::ModelOptions::option (util::MAX_DENS_CORRECTION))
	infStepMaxDens = std::max(infStepMaxDens, timeStepMaxDensity);
    timeStepMaxDensity = infStepMaxDens;
    
    totalDensity += (*iter)->getDensity();
    if ((*iter)->getStartDate() == (Global::simulationTime-1)) {
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


// -----  Data checkpointing  -----

void DescriptiveWithinHostModel::checkpoint (istream& stream) {
    WithinHostModel::checkpoint (stream);
    for(int i=0; i<_MOI; ++i) {
	loadInfection(stream);	// create infections using a virtual function call
    }
}
void DescriptiveWithinHostModel::checkpoint (ostream& stream) {
    WithinHostModel::checkpoint (stream);
    BOOST_FOREACH (DescriptiveInfection* inf, infections) {
	(*inf) & stream;
    }
}

} }