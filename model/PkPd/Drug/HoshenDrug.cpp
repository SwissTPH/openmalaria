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

#include "PkPd/Drug/HoshenDrug.h"

#include <cmath>
#include <assert.h>

namespace OM { namespace PkPd {
    
// -----  Static variables and functions  -----

double HoshenDrug::minutesPerTimeStep;

void HoshenDrug::init () {
  minutesPerTimeStep = Global::interval * 24*60;
}


// -----  Non-static variables and functions  -----

HoshenDrug::HoshenDrug(const HoshenDrugType* type) :
    typeData (type),
    _concentration (0),
    _nextConcentration (0)
{}


void HoshenDrug::addDose (double concentration, int delay) {
    assert (delay == 0);
    _concentration += concentration;
    _nextConcentration = _concentration * decayFactor (minutesPerTimeStep);
}

bool HoshenDrug::decay() {
    _concentration = _nextConcentration;
    _nextConcentration = _concentration * decayFactor (minutesPerTimeStep);
    //TODO: if concentration is negligible, return true to clean up this object
    return false;
}

double HoshenDrug::calculateDrugFactor(uint32_t proteome_ID) const {
  //Returning an average of 2 points
  double param = ((HoshenDrugType*)typeData)->proteomePDParameters.find(proteome_ID)->second;
  double startFactor = 3.8/(1+param/_concentration);
  double endFactor = 3.8/(1+param/_nextConcentration);
  return std::exp(-(startFactor + endFactor)/2);
}

double HoshenDrug::decayFactor (double time) {
  //k = log(2)/halfLife
  return std::exp(-time*log(2.0)/((HoshenDrugType*)typeData)->halfLife);
}

} }