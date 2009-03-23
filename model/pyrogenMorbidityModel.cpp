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

#include "pyrogenMorbidityModel.h"
#include "inputData.h"

using namespace std;

double PyrogenMorbidityModel::initPyroThres;
double PyrogenMorbidityModel::smuY;
double PyrogenMorbidityModel::Ystar2_13;
double PyrogenMorbidityModel::alpha14;
double PyrogenMorbidityModel::Ystar1_26;

void PyrogenMorbidityModel::init(){
  initPyroThres=getParameter(Params::Y_STAR_0);
  smuY=-log(0.5)/(daysInYear/Global::interval*getParameter(Params::Y_STAR_HALF_LIFE));
  Ystar2_13=getParameter(Params::Y_STAR_SQ);
  alpha14=getParameter(Params::ALPHA);
  Ystar1_26=getParameter(Params::Y_STAR_1);
}

PyrogenMorbidityModel::PyrogenMorbidityModel(double cF) :
     MorbidityModel (cF), _pyrogenThres (initPyroThres)
{}

double PyrogenMorbidityModel::getPEpisode(double timeStepMaxDensity, double totalDensity) {
  updatePyrogenThres(totalDensity);
  return 1-1/(1+(timeStepMaxDensity/_pyrogenThres));;
}

double PyrogenMorbidityModel::getPyrogenThres(){
  return _pyrogenThres;
}

void PyrogenMorbidityModel::updatePyrogenThres(double totalDensity){
  int i;
  //Number of categories in the numerical approx. below
  const int n= 11;
  const double delt= 1.0/n;
  double valYstar=_pyrogenThres;
  //Numerical approximation to equation 2, AJTMH p.57
  for ( i=1;i<=n; i++) {
    valYstar=valYstar+totalDensity*alpha14*Global::interval*delt/((Ystar1_26+totalDensity)*(Ystar2_13+valYstar))-smuY*valYstar*delt;
  }
  _pyrogenThres=valYstar;
}

void PyrogenMorbidityModel::read(istream& in) {
  in >> _comorbidityFactor;
  in >> _pyrogenThres;
}

void PyrogenMorbidityModel::write(ostream& out) const {
  out << _comorbidityFactor << endl;
  out << _pyrogenThres << endl;
}
