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

#include "WithinHost/Pathogenesis/Submodels.h"
#include "Parameters.h"

#include <cmath>
using namespace std;

namespace OM {
namespace WithinHost {
namespace Pathogenesis {
    using namespace Monitoring;

// ———  Müller presentation model  ———

double MuellerPathogenesis::rateMultiplier_31;
double MuellerPathogenesis::densityExponent_32;

void MuellerPathogenesis::init( const Parameters& parameters ){
  rateMultiplier_31=parameters[Parameters::MUELLER_RATE_MULTIPLIER];
  densityExponent_32=parameters[Parameters::MUELLER_DENSITY_EXPONENT];
}

double MuellerPathogenesis::getPEpisode(double, double totalDensity) {
  double incidenceDensity = rateMultiplier_31 * (pow(totalDensity, densityExponent_32)) * TimeStep::yearsPerInterval;
  return 1.0-exp(-incidenceDensity);
}


// ———  Pyrogenic threshold model  ———

double PyrogenPathogenesis::initPyroThres;
double PyrogenPathogenesis::smuY;
double PyrogenPathogenesis::Ystar2_13;
double PyrogenPathogenesis::alpha14;
double PyrogenPathogenesis::Ystar1_26;

void PyrogenPathogenesis::init( const Parameters& parameters ){
  initPyroThres=parameters[Parameters::Y_STAR_0];
  smuY=-log(0.5)/(TimeStep::stepsPerYear*parameters[Parameters::Y_STAR_HALF_LIFE]);
  Ystar2_13=parameters[Parameters::Y_STAR_SQ];
  alpha14=parameters[Parameters::ALPHA];
  Ystar1_26=parameters[Parameters::Y_STAR_1];
}

PyrogenPathogenesis::PyrogenPathogenesis(double cF) :
     PathogenesisModel (cF), _pyrogenThres (initPyroThres)
{}


double PyrogenPathogenesis::getPEpisode(double timeStepMaxDensity, double totalDensity) {
  updatePyrogenThres(totalDensity);
  return 1-1/(1+(timeStepMaxDensity/_pyrogenThres));;
}

void PyrogenPathogenesis::summarize (const Host::Human& human) {
    Survey::current()
        .addDouble( Report::MD_PYROGENIC_THRESHOLD, human, _pyrogenThres )
        .addDouble( Report::MD_LOG_PYROGENIC_THRESHOLD, human, log(_pyrogenThres+1.0) );
}

void PyrogenPathogenesis::updatePyrogenThres(double totalDensity){
    // Note: this calculation is slow (something like 5% of runtime)
    
  //Number of categories in the numerical approx. below
  const int n= 11;
  const double delt= 1.0/n;
  //Numerical approximation to equation 2, AJTMH p.57
  for (int i=1;i<=n; ++i) {
    _pyrogenThres += totalDensity * alpha14 * TimeStep::interval * delt /
        ( (Ystar1_26 + totalDensity) * (Ystar2_13 + _pyrogenThres) )
        - smuY * _pyrogenThres * delt;
  }
}


void PyrogenPathogenesis::checkpoint (istream& stream) {
    PathogenesisModel::checkpoint (stream);
    _pyrogenThres & stream;
}
void PyrogenPathogenesis::checkpoint (ostream& stream) {
    PathogenesisModel::checkpoint (stream);
    _pyrogenThres & stream;
}


// ———  Predetermined episodes model  ———

double PredetPathogenesis::getPEpisode(double timeStepMaxDensity, double totalDensity) {
  updatePyrogenThres(totalDensity);
  if ( timeStepMaxDensity > _pyrogenThres) {
      return 1;
    }
  else{
    return 0;
  }
}

} } }
