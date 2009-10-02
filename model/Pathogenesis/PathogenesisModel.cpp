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

#include "Pathogenesis/PathogenesisModel.h"
#include "Pathogenesis/Pyrogen.h"
#include "Pathogenesis/Predet.h"
#include "Pathogenesis/Mueller.h"
#include "WithinHost/WithinHostModel.h"
#include "inputData.h"
#include "util/gsl.h"

using namespace std;


//BEGIN static
double PathogenesisModel::indirRiskCoFactor_18;
double PathogenesisModel::sevMal_21;
double PathogenesisModel::comorbintercept_24;
double PathogenesisModel::critAgeComorb_30;


void PathogenesisModel::init() {
  indirRiskCoFactor_18=(1-exp(-getParameter(Params::INDIRECT_RISK_COFACTOR)));
  sevMal_21=getParameter(Params::SEVERE_MALARIA_THRESHHOLD);
  comorbintercept_24=1-exp(-getParameter(Params::COMORBIDITY_INTERCEPT));
  critAgeComorb_30=getParameter(Params::CRITICAL_AGE_FOR_COMORBIDITY);
  
  if (Global::modelVersion & PREDETERMINED_EPISODES) {
    //no separate init:
    PyrogenPathogenesis::init();
  } else {
    if (Global::modelVersion & MUELLER_PRESENTATION_MODEL)
      MuellerPathogenesis::init();
    else
      PyrogenPathogenesis::init();
  }
}

PathogenesisModel* PathogenesisModel::createPathogenesisModel(double cF) {
  if (Global::modelVersion & PREDETERMINED_EPISODES) {
    return new PredetPathogenesis(cF);
  }
  else {
    if (Global::modelVersion & MUELLER_PRESENTATION_MODEL) {
      return new MuellerPathogenesis(cF);
    }
    else {
      return new PyrogenPathogenesis(cF);
    }
  }
}

PathogenesisModel* PathogenesisModel::createPathogenesisModel(istream& in) {
  if (Global::modelVersion & PREDETERMINED_EPISODES) {
    return new PredetPathogenesis(in);
  }
  else {
    if (Global::modelVersion & MUELLER_PRESENTATION_MODEL) {
      return new MuellerPathogenesis(in);
    }
    else {
      return new PyrogenPathogenesis(in);
    }
  }
}
//END static


PathogenesisModel::PathogenesisModel(double cF) :
    _comorbidityFactor(cF)
{}
PathogenesisModel::PathogenesisModel(istream& in)
{
  in >> _comorbidityFactor;
}

Pathogenesis::State PathogenesisModel::determineState (double ageYears, WithinHostModel& withinHostModel) {
  double timeStepMaxDensity = withinHostModel.getTimeStepMaxDensity();
  double prEpisode = getPEpisode(timeStepMaxDensity, withinHostModel.getTotalDensity());
  
  //Decide whether a clinical episode occurs and if so, which type
  double pCoinfection=comorbintercept_24/(1+ageYears/critAgeComorb_30);
  pCoinfection*=_comorbidityFactor;
  
  if ((gsl::rngUniform()) < prEpisode) {
    //Fixed severe threshold
    double severeMalThreshold=sevMal_21+1;
    double prSevereEpisode=1-1/(1+timeStepMaxDensity/severeMalThreshold);
    
    Pathogenesis::State ret;
    if (gsl::rngUniform() < prSevereEpisode)
      ret = Pathogenesis::STATE_SEVERE;
    else if (gsl::rngUniform() < pCoinfection)
      ret = Pathogenesis::STATE_COINFECTION;
    else
      ret = Pathogenesis::STATE_MALARIA;
    
    /* Indirect mortality	
       IndirectRisk is the probability of dying from indirect effects of malaria
       conditional on not having an acute attack of malaria
    */
    double indirectRisk=indirRiskCoFactor_18/(1+ageYears/critAgeComorb_30);
    indirectRisk*=_comorbidityFactor;
    if (gsl::rngUniform() < indirectRisk)
      ret = Pathogenesis::State (ret | Pathogenesis::INDIRECT_MORTALITY);
    
    return ret;
  } else if(Global::modelVersion & NON_MALARIA_FEVERS) {
    //TODO: should this be stored in the XML file?
    const double RelativeRiskNonMalariaFever= 1.0;
    double prNonMalariaFever=pCoinfection*RelativeRiskNonMalariaFever;
    if ((gsl::rngUniform()) < prNonMalariaFever)
      return Pathogenesis::SICK;
  }
  return Pathogenesis::NONE;
}

void PathogenesisModel::write(ostream& out) const {
  out << _comorbidityFactor << endl;
}
