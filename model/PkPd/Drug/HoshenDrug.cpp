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

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

HoshenDrug::HoshenDrug(const DrugType* type) : Drug(type) {
}

HoshenDrug::HoshenDrug (const DrugType* type, istream& in) :
  Drug(type, in) {
}

void HoshenDrug::write (ostream& out) const {
    Drug::write(out);
}



double HoshenDrug::calculateDrugFactor(const ProteomeInstance* infProteome) const {
  //Returning an average of 2 points
  double param = typeData->proteomePDParameters.find(infProteome->getProteomeID())->second;
  double startFactor = 3.8/(1+param/_concentration);
  double endFactor = 3.8/(1+param/_nextConcentration);
  return exp(-(startFactor + endFactor)/2);
}

double HoshenDrug::decayFactor (double time) {
  //k = log(2)/halfLife
  return exp(-time*log(2.0)/typeData->halfLife);
}

