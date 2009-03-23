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

#include "morbidityModel.h"
#include "pyrogenMorbidityModel.h"
#include "predetMorbidityModel.h"
#include "muellerMorbidityModel.h"
#include "global.h"
#include "inputData.h"
#include "GSLWrapper.h"

using namespace std;


double MorbidityModel::indirRiskCoFactor_18;
double MorbidityModel::sevMal_21;
double MorbidityModel::comorbintercept_24;
double MorbidityModel::critAgeComorb_30;


void MorbidityModel::initModels() {
  indirRiskCoFactor_18=(1-exp(-getParameter(Params::INDIRECT_RISK_COFACTOR)));
  sevMal_21=getParameter(Params::SEVERE_MALARIA_THRESHHOLD);
  comorbintercept_24=1-exp(-getParameter(Params::COMORBIDITY_INTERCEPT));
  critAgeComorb_30=getParameter(Params::CRITICAL_AGE_FOR_COMORBIDITY);
  
  
  // Perhaps not all models are used, but currently this initialization is fast
  // anyway:
  
  PyrogenMorbidityModel::init();
  //no separate init:
  //PredetMorbidityModel::init();
  MuellerMorbidityModel::init();
}

MorbidityModel* MorbidityModel::createMorbidityModel(double cF) {
  if (Global::modelVersion & PREDETERMINED_EPISODES) {
    return new PredetMorbidityModel(cF);
  }
  else {
    if (Global::modelVersion & MUELLER_MORBIDITY_MODEL) {
      return new MuellerMorbidityModel(cF);
    }
    else {
      return new PyrogenMorbidityModel(cF);
    }
  }
}


MorbidityModel::MorbidityModel(double cF) :
    _comorbidityFactor(cF)
{}

Morbidity::Infection MorbidityModel::infectionEvent(double ageYears, double totalDensity, double timeStepMaxDensity) {
  double prEpisode = getPEpisode(timeStepMaxDensity,totalDensity);
  
  //Decide whether a clinical episode occurs and if so, which type
  double pCoinfection=comorbintercept_24/(1+ageYears/critAgeComorb_30);
  pCoinfection*=_comorbidityFactor;
  
  if ((W_UNIFORM()) < prEpisode) {
    //Fixed severe threshold
    double severeMalThreshold=sevMal_21+1;
    double prSevereEpisode=1-1/(1+timeStepMaxDensity/severeMalThreshold);
    
    Morbidity::Infection ret = Morbidity::UNCOMPLICATED;
    
    if (W_UNIFORM() < prSevereEpisode)
      ret = Morbidity::SEVERE;
    else if (W_UNIFORM() < pCoinfection)
      ret = Morbidity::COINFECTION;
    
    /*
      Indirect mortality	
      IndirectRisk is the probability of dying from indirect effects of malaria
      conditional on not having an acute attack of malaria
    *
    double indirectRisk=indirRiskCoFactor_18/(1+ageYears/critAgeComorb_30);
    indirectRisk*=_comorbidityFactor;
    if (W_UNIFORM() < indirectRisk)
      ret = Morbidity::Infection (ret | Morbidity::INDIRECT_MORTALITY);
    */
    return ret;
  }
  else if(Global::modelVersion & NON_MALARIA_FEVERS) {
    double prNonMalariaFever=pCoinfection*RelativeRiskNonMalariaFever;
    if ((W_UNIFORM()) < prNonMalariaFever)
      return Morbidity::NON_MALARIA;
  }
  return Morbidity::NONE;
}

bool MorbidityModel::indirectDeath(double ageYears) {
  double indirectRisk=indirRiskCoFactor_18/(1+ageYears/critAgeComorb_30);
  indirectRisk*=_comorbidityFactor;
  return (W_UNIFORM() < indirectRisk);
}

double MorbidityModel::getPyrogenThres(){
  return 0;
}


void MorbidityModel::read(istream& in) {
  in >> _comorbidityFactor;
}

void MorbidityModel::write(ostream& out) const {
  out << _comorbidityFactor << endl;
}
