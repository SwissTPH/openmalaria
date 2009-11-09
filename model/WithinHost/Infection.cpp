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

float Infection::cumulativeYstar;
float Infection::cumulativeHstar;
double Infection::alpha_m;
double Infection::decayM;

void Infection::init () {
  cumulativeYstar=(float)getParameter(Params::CUMULATIVE_Y_STAR);
  cumulativeHstar=(float)getParameter(Params::CUMULATIVE_H_STAR);
  alpha_m=1-exp(-getParameter(Params::NEG_LOG_ONE_MINUS_ALPHA_M));
  decayM=getParameter(Params::DECAY_M);
}


Infection::Infection (istream& in) {
  in >> _startdate;
  in >> _density;
  if (Global::modelVersion & INCLUDES_PK_PD) {
    int proteomeID;
    in >> proteomeID;
    _proteome = ProteomeInstance::getProteome(proteomeID);
  } else
    _proteome = NULL;
}
void Infection::write (ostream& out) const {
  out << _startdate << endl; 
  out << _density << endl; 
  if (Global::modelVersion & INCLUDES_PK_PD) {
    out << _proteome->getProteomeID() << endl; 
  }
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
