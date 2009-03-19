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

#include "muellerMorbidityModel.h"
#include "inputData.h"
using namespace std;

double MuellerMorbidityModel::rateMultiplier_31;
double MuellerMorbidityModel::densityExponent_32;

void MuellerMorbidityModel::init(){
  rateMultiplier_31=getParameter(Params::MUELLER_RATE_MULTIPLIER);
  densityExponent_32=getParameter(Params::MUELLER_DENSITY_EXPONENT);
}

double MuellerMorbidityModel::getPEpisode(double timeStepMaxDensity, double totalDensity) {
 double incidenceDensity=rateMultiplier_31*(pow(totalDensity, densityExponent_32))/(1.0*Global::intervalsPerYear);
   return 1-exp(-incidenceDensity);
  return 0;
}

void MuellerMorbidityModel::read(istream& in) {
  //Empty
}

void MuellerMorbidityModel::write(ostream& out) const {
 //Empty
}
