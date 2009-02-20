/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include "caseManagement.h"
#include "inputData.h"
#include "global.h"

CaseManagement::CaseManagement(){

  _oddsRatioThreshold=exp(get_parameter(17));
  for (int regimen=0;regimen<3; regimen++) {
    _cureRate[regimen]=get_curerate(regimen);
    _probGetsTreatment[regimen]=get_p_gets_treatment(regimen);
    _probParasitesCleared[regimen]=get_p_parasites_cleared(regimen);
  }
  for (int agegrp=0;agegrp<2; agegrp++) {
    _probSequelaeTreated[agegrp]=get_p_sequelae(agegrp);
    _probSequelaeUntreated[agegrp]=get_p_sequelae(agegrp);
  }
  _caseManagementMemory=get_health_system_memory();
  readCaseFatalityRatio();
}

CaseManagement::~CaseManagement(){
}

void CaseManagement::readCaseFatalityRatio(){
  int numOfGroups=get_number_of_cfrgroups();
  _inputAge.resize(numOfGroups+1);
  _caseFatalityRate.resize(numOfGroups+1);
  for (int i=0;i<numOfGroups; i++) {
    _inputAge[i]=get_cfr_lb(i);
    _caseFatalityRate[i]=get_cfr(i);
  }
  //Use a high number for the upper bound.
  _inputAge[numOfGroups]=999;
  //CFR is constant for everyone above the highest upperbound
  _caseFatalityRate[numOfGroups]=get_cfr(numOfGroups-1);
  if (( numOfGroups == 1) && ( _caseFatalityRate[0] ==  0)) {
    _noMortality=true;
  }
  else {
    _noMortality=false;
  }
}

double CaseManagement::getCommunityCaseFatalityRate(double caseFatalityRatio) const{
  return caseFatalityRatio*_oddsRatioThreshold/(1-caseFatalityRatio+caseFatalityRatio*_oddsRatioThreshold);
    
}

int CaseManagement::getNextRegimen(int simulationTime, int diagnosis, int tLastTreated, int regimen) const{
  if (diagnosis == Diagnosis::SEVERE_MALARIA){
    return 3;
  }
  else if (tLastTreated > ( simulationTime-_caseManagementMemory)){
    return 2;
  }
  return 1;
}

int CaseManagement::getCaseManagementMemory() const{
  return _caseManagementMemory;
}

double CaseManagement::getRiskFromMaternalInfection() const{
  return _riskFromMaternalInfection;
}

void CaseManagement::setRiskFromMaternalInfection(double risk){
  _riskFromMaternalInfection=risk;
}

double CaseManagement::getProbabilityGetsTreatment(int regimen) const{
  return _probGetsTreatment[regimen];
}

double CaseManagement::getProbabilityParasitesCleared(int regimen) const{
  return _probParasitesCleared[regimen];
}

double CaseManagement::getCureRate(int regimen) const{
  return _cureRate[regimen];
}

double CaseManagement::getProbabilitySequelaeTreated(int ageGroup) const{
  return _probSequelaeTreated[ageGroup];
}

double CaseManagement::getProbabilitySequelaeUntreated(int ageGroup) const{
  return _probSequelaeUntreated[ageGroup];
}

double CaseManagement::caseFatality(double ageYears){
    /*
    Linear interpolation to get age-specific hospital case fatality rates
    ageyears: age of person in years
    */

    double a0=0.0;
    double a1=0.0;
    double f0=0.0;
    double f1=0.0;
    int i;
    i=1;
    if (_noMortality) {
        return 0;
    }
    else {
        while( _inputAge[i - 1] <=  ageYears) {
            i=i+1;
            a0=_inputAge[i-1 - 1];
            a1=_inputAge[i - 1];
            f0=_caseFatalityRate[i-1 - 1];
            f1=_caseFatalityRate[i - 1];
        }
        return (ageYears-a0)/(a1-a0)*(f1-f0)+f0;
    }
}

