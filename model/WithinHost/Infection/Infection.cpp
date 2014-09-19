/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "WithinHost/Infection/Infection.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"

#include <cmath>

namespace OM { namespace WithinHost {
    
double Infection::invCumulativeYstar;
double Infection::invCumulativeHstar;
double Infection::alpha_m;
double Infection::decayM;
SimTime Infection::latentP;

void Infection::init (const OM::Parameters& parameters, SimTime latP) {
    latentP = latP;
    // calculate inverses here, so we can use multiplication later (faster):
    invCumulativeYstar = 1.0 / parameters[Parameters::CUMULATIVE_Y_STAR];
    invCumulativeHstar = 1.0 / parameters[Parameters::CUMULATIVE_H_STAR];
    alpha_m = 1.0 - exp(-parameters[Parameters::NEG_LOG_ONE_MINUS_ALPHA_M]);
    decayM = parameters[Parameters::DECAY_M];
}


double Infection::immunitySurvivalFactor (double ageInYears, double cumulativeh, double cumulativeY) {
  //Documentation: AJTMH pp22-23
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
    dH=1.0 / (1.0 + (cumulativeh-1.0) * invCumulativeHstar);
    dY=1.0 / (1.0 + (cumulativeY - m_cumulativeExposureJ) * invCumulativeYstar);
  }
  dA = 1.0 - alpha_m * exp(-decayM * ageInYears);
  double ret = std::min(dY*dH*dA, 1.0);
  util::streamValidate( ret );
  return ret;
}


Infection::Infection (istream& stream) : m_startDate(sim::never()) {
    m_startDate & stream;
    m_proteome_ID & stream;
    m_density & stream;
    m_cumulativeExposureJ & stream;
}
void Infection::checkpoint (ostream& stream) {
    m_startDate & stream;
    m_proteome_ID & stream;
    m_density & stream;
    m_cumulativeExposureJ & stream;
}

} }