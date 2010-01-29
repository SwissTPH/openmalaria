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

#include "WithinHost/Infection.h"
#include "inputData.h"
#include "util/ModelOptions.hpp"

#include <cmath>

namespace OM { namespace WithinHost {
    
float Infection::cumulativeYstar;
float Infection::cumulativeHstar;
double Infection::alpha_m;
double Infection::decayM;

void Infection::init () {
  cumulativeYstar=(float)InputData.getParameter(Params::CUMULATIVE_Y_STAR);
  cumulativeHstar=(float)InputData.getParameter(Params::CUMULATIVE_H_STAR);
  alpha_m=1-exp(-InputData.getParameter(Params::NEG_LOG_ONE_MINUS_ALPHA_M));
  decayM=InputData.getParameter(Params::DECAY_M);
}


double Infection::immunitySurvivalFactor (double ageInYears, double cumulativeh, double cumulativeY) {
  //effect of cumulative Parasite density (named Dy in AJTM)
  double dY;
  //effect of number of infections experienced since birth (named Dh in AJTM)
  double dH;
  //effect of age-dependent maternal immunity (named Dm in AJTM)
  double dA;
  
  if (cumulativeh <= 1.0) {
    dY=1.0;
    dH=1.0;
  } else {
    dH=1.0 / (1.0 + (cumulativeh-1.0) / cumulativeHstar);
    //TODO: compare this with the asex paper
    dY=1.0 / (1.0 + (cumulativeY-_cumulativeExposureJ) / cumulativeYstar);
  }
  dA = 1.0 - alpha_m * exp(-decayM * ageInYears);
  return std::min(dY*dH*dA, 1.0);
}


void Infection::checkpoint (istream& stream) {
    _startdate & stream;
    _density & stream;
    _cumulativeExposureJ & stream; 
    proteome_ID & stream;
}
void Infection::checkpoint (ostream& stream) {
    _startdate & stream;
    _density & stream;
    _cumulativeExposureJ & stream; 
    proteome_ID & stream;
}

} }