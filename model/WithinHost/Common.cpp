/*
  This file is part of OpenMalaria.
 
  Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "WithinHost/Common.h"
#include "inputData.h"
#include "util/errors.hpp"

using namespace std;

namespace OM { namespace WithinHost {

CommonInfection* (* CommonWithinHost::createInfection) (uint32_t protID);
CommonInfection* (* CommonWithinHost::checkpointedInfection) (istream& stream);


// -----  Initialization  -----

CommonWithinHost::CommonWithinHost() :
    WithinHostModel(), pkpdModel(PkPd::PkPdModel::createPkPdModel ())
{}
CommonWithinHost::~CommonWithinHost() {
    delete pkpdModel;
    for(std::list<CommonInfection*>::iterator it = infections.begin(); it != infections.end(); ++it) {
	delete *it;
    }
}


// -----  Simple infection adders/removers  -----

void CommonWithinHost::newInfection(){
    if (_MOI < MAX_INFECTIONS) {
	infections.push_back(createInfection (pkpdModel->new_proteome_ID ()));
	_MOI++;
	_cumulativeh++;
    }
}

void CommonWithinHost::clearAllInfections(){
    for(std::list<CommonInfection*>::iterator it = infections.begin(); it != infections.end(); ++it) {
	delete *it;
    }
    infections.clear();
    _MOI=0;
}


// -----  medicate drugs -----

void CommonWithinHost::medicate(string drugName, double qty, double time, double age) {
  pkpdModel->medicate(drugName, qty, time, age);
}


// -----  Density calculations  -----

void CommonWithinHost::calculateDensities(double ageInYears, double BSVEfficacy) {
    updateImmuneStatus ();
    
    totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  double timestepCumY = _cumulativeY;
  
  for(std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end();){
    double survivalFactor = (1.0-BSVEfficacy) * _innateImmSurvFact;
    survivalFactor *= pkpdModel->getDrugFactor((*inf)->get_proteome_ID());
    survivalFactor *= (*inf)->immunitySurvivalFactor(ageInYears, _cumulativeh, timestepCumY);
    
    // We update the density, and if updateDensity returns true (parasites extinct) then remove the infection.
    if ((*inf)->updateDensity(Global::simulationTime, survivalFactor)) {
	delete *inf;
      inf = infections.erase(inf);	// inf points to next infection now so don't increment with ++inf
      --_MOI;
      continue;	// infection no longer exists so skip the rest
    }
    
    totalDensity += (*inf)->getDensity();
    timeStepMaxDensity = max(timeStepMaxDensity, (*inf)->getDensity());
    _cumulativeY += Global::interval*(*inf)->getDensity();
    
    ++inf;
  }
  
  pkpdModel->decayDrugs ();
}


// -----  Summarize  -----

int CommonWithinHost::countInfections (int& patentInfections) {
    if (infections.empty()) return 0;
    patentInfections = 0;
    for (std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end(); ++inf) {
	if ((*inf)->getDensity() > detectionLimit)
	    patentInfections++;
    }
    return infections.size();
}


void CommonWithinHost::checkpoint (istream& stream) {
    WithinHostModel::checkpoint (stream);
    (*pkpdModel) & stream;
    for (int i = 0; i < _MOI; ++i) {
	infections.push_back (checkpointedInfection (stream));
    }
}
void CommonWithinHost::checkpoint (ostream& stream) {
    WithinHostModel::checkpoint (stream);
    (*pkpdModel) & stream;
    for(std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end();) {
	(**inf) & stream;
    }
}

} }