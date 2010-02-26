/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "Host/NeonatalMortality.h"
#include "Host/Human.h"
#include "util/random.h"

namespace OM { namespace Host {
    using namespace OM::util;

double NeonatalMortality::_riskFromMaternalInfection = 0.0;
std::vector<double> NeonatalMortality::_prevalenceByGestationalAge;


void NeonatalMortality::init() {
  int timeStepsPer5Months = 150 / Global::interval;
  _prevalenceByGestationalAge.assign(timeStepsPer5Months, 0.0);
}

void NeonatalMortality::staticCheckpoint (istream& stream) {
    _riskFromMaternalInfection & stream;
    _prevalenceByGestationalAge & stream;
}
void NeonatalMortality::staticCheckpoint (ostream& stream) {
    _riskFromMaternalInfection & stream;
    _prevalenceByGestationalAge & stream;
}

bool NeonatalMortality::eventNeonatalMortality() {
  return random::uniform_01() <= _riskFromMaternalInfection;
}

void NeonatalMortality::update (const list<Host::Human>& population) {
  // For individuals in the age range 20-25, we sum:
  int nCounter=0;	// total number
  int pCounter=0;	// number with patent infections, needed for prev in 20-25y
  
  for (std::list<Host::Human>::const_iterator iter = population.begin(); iter != population.end(); ++iter){
    double ageYears = iter->getAgeInYears();
    // Note: since we're using a linked list, we have to iterate until we reach
    // the individuals we're interested in. Due to population structure, it's
    // probably quickest to start iterating from the oldest.
    if(ageYears >= 25.0) continue;
    if (ageYears < 20.0) break;	// Not interested in younger individuals.
    
    nCounter ++;
    if (iter->detectibleInfection())
      pCounter ++;
  }
  
  NeonatalMortality::calculateRiskFromMaternalInfection(nCounter, pCounter);
}

void NeonatalMortality::calculateRiskFromMaternalInfection (int nCounter, int pCounter) {
  //Goodman estimated for neonatal mortality due to malaria in pregnancy
  const double gEst = 0.011;
  //Critical value of Prev20-25 for neonatal mortality
  const double critPrev2025 = 0.25;
  //Critical value for estimating prevalence in primigravidae
  const double critPrevPrim = 0.19;
  //Proportion of births with primigravid mothers
  const double pBirthPrim = 0.3;
  //default value for prev2025, for short simulations 
  double prev2025 = 0.25;
  prev2025 = double(pCounter) / nCounter;  
  double maxprev = prev2025;
  //gestational age is in time steps for the last 5 months of pregnancy only
  int timeStepsMinus1 = 150 / Global::interval - 1;
  //update the vector containing the prevalence by gestational age
  for (int t=0; t < timeStepsMinus1; t++) {
    _prevalenceByGestationalAge[t] = _prevalenceByGestationalAge[t+1];
    if (_prevalenceByGestationalAge[t] > maxprev) {
      maxprev = _prevalenceByGestationalAge[t];
    }
  }
  _prevalenceByGestationalAge[timeStepsMinus1] = prev2025;
  //equation (2) p 75 AJTMH 75 suppl 2
  double prevpg= maxprev / (critPrevPrim + maxprev);
  //equation (1) p 75 AJTMH 75 suppl 2
  _riskFromMaternalInfection = gEst * pBirthPrim * (1.0-exp(-prevpg/critPrev2025));
}

} }