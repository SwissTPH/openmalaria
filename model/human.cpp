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
#include <string>
#include <string.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include "human.h"
#include "simulation.h"
#include "inputData.h"
#include "GSLWrapper.h"
#include "caseManagement.h"
#include "summary.h"
#include "proteome.h"
#include "intervention.h"
#include "transmissionModel.h"
#include "oldWithinHostModel.h"


/*
  Constants common to all humans
  Decay in anti-toxic immunity
*/
//TODOConversion
#ifdef _WIN32
#define isnan(x) ((x) != (x))
#endif

double smuY;
double sigma_i;
double initPyroThres;
double Ystar2_13;
double alpha14;
double Ystar1_26;
double sevMal_21;
double critAgeComorb_30;
double comorbintercept_24;
double indirRiskCoFactor_18;
double immPenalty_22;
double asexImmRemain;
double immEffectorRemain;
double rateMultiplier_31;
double densityExponent_32;
double BaselineAvailabilityShapeParam;
double detectionlimit;

void initHumanParameters () {
  /*

    Init parameters that are common to all humans

  */

  double densitybias;
  smuY=-log(0.5)/(daysInYear/interval*get_parameter(25));
  sigma_i=sqrt(get_parameter(6));
  initPyroThres=get_parameter(28);
  Ystar2_13=get_parameter(13);
  alpha14=get_parameter(14);
  indirRiskCoFactor_18=(1-exp(-get_parameter(18)));
  sevMal_21=get_parameter(21);
  immPenalty_22=1-exp(get_parameter(22));
  comorbintercept_24=1-exp(-get_parameter(24));
  Ystar1_26=get_parameter(26);
  immEffectorRemain=exp(-get_parameter(23));
  asexImmRemain=exp(-get_parameter(27));
  critAgeComorb_30=get_parameter(30);
  if ( isOptionIncluded(modelVersion, MuellerMorbidityModel)) {
    rateMultiplier_31=get_parameter(31);
    densityExponent_32=get_parameter(32);
  }
  /*
    TODO: This densitiybias function should be part of the scenario description XML, not the parameter element.
    or maybe it should be a parameter, as we want to fit it... but the garki analysis numbers are a bit dangerous
    add an attribute to scenario.xml densityQuantification="malariaTherapy|garki|other"
  */
  if (( get_analysis_no() <  22) || ( get_analysis_no() >  30)) {
    densitybias=get_parameter(15);
  }
  else {
    densitybias=get_parameter(20);
  }
  detectionlimit=get_detectionlimit()*densitybias;
  BaselineAvailabilityShapeParam=get_parameter(16);
}

Human::Human() {
  _withinHostModel = new OldWithinHostModel(this);
  if (isOptionIncluded(modelVersion, includesPKPD)) {
    _drugs = list<Drug*>();
    _proxy = new DrugProxy(this);
  }
}


Human::Human(int ID, int dateOfBirth, CaseManagement* caseManagement, int simulationTime) 
  : _simulationTime(simulationTime), _caseManagement(caseManagement) {
 
  int i;
  //std::cout<<"newH:ID dateOfBirth "<<ID<<" "<<dateOfBirth<<std::endl;
  _BSVEfficacy=0.0;
  _cumulativeEIRa=0.0;
  _cumulativeh=0.0;
  _cumulativeInfections=0;
  _cumulativeY=0.0;
  _cumulativeYlag=0.0;
  _dateOfBirth=dateOfBirth;
  _doomed=0;
  _timeStepMaxDensity=0.0;
  _ID=ID;
  _innateImmunity=(double)(W_GAUSS((0), (sigma_i)));
  _lastVaccineDose=0;
  _latestRegimen=0;
  _tLastTreatment=missing_value;
  _MOI=0;
  _patentinfections=0;
  _PEVEfficacy=0.0;
  _pinfected=0.0;
  _ptransmit=0.0;
  _pyrogenThres=initPyroThres;
  _TBVEfficacy=0.0;
  _totalDensity=0.0;
  for ( i=1;i<=4; i++) {
    _ylag[i-1]=0.0;
  }
  _ITN=false;/* bool */;
  _latestEvent.setTime(missing_value);
  _latestEvent.setCaseManagement(caseManagement);
  _lastSPDose=missing_value;
  _lastIptiOrPlacebo=missing_value;
  if (isOptionIncluded(modelVersion, negativeBinomialMassAction)) {
    _BaselineAvailabilityToMosquitoes=(W_GAMMA((BaselineAvailabilityShapeParam), (BaselineAvailabilityMean/BaselineAvailabilityShapeParam)));
  }
  else if( isOptionIncluded(modelVersion, lognormalMassAction)) {
    _BaselineAvailabilityToMosquitoes=(W_LOGNORMAL((log(BaselineAvailabilityMean))-(0.5*pow(BaselineAvailabilityShapeParam, 2)), (BaselineAvailabilityShapeParam)));
  }
  else if (isOptionIncluded(modelVersion, transHet)) {
    _BaselineAvailabilityToMosquitoes=0.2;
    if (W_UNIFORM() < 0.5) {            
      _BaselineAvailabilityToMosquitoes=1.8;
    }
  }
  else {
    _BaselineAvailabilityToMosquitoes=BaselineAvailabilityMean;
  }
  if (isOptionIncluded(modelVersion, comorbHet)) {
    _comorbidityFactor=0.2;
    if (W_UNIFORM() < 0.5) {            
      _comorbidityFactor=1.8;
    }	
  }
  else {
    _comorbidityFactor=1.0;
  }
  if (isOptionIncluded(modelVersion, treatHet)) {
    _treatmentSeekingFactor=0.2;
    if (W_UNIFORM() < 0.5) {            
      _treatmentSeekingFactor=1.8;
    }	
  }
  else {
    _treatmentSeekingFactor=1.0;
  }
  if (isOptionIncluded(modelVersion, transTreatHet)) {
    _treatmentSeekingFactor=0.2;
    _BaselineAvailabilityToMosquitoes=1.8;
    if (W_UNIFORM()<0.5) {
      _treatmentSeekingFactor=1.8;
      _BaselineAvailabilityToMosquitoes=0.2;
    }
  }
  if (isOptionIncluded(modelVersion, comorbTreatHet)) {
    _comorbidityFactor=0.2;
    _treatmentSeekingFactor=1.8;
    if (W_UNIFORM()<0.5) {
      _treatmentSeekingFactor=0.2;
      _comorbidityFactor=1.8;
    }
  }
  if (isOptionIncluded(modelVersion,comorbTransHet)) {
    _BaselineAvailabilityToMosquitoes=1.8;
    _comorbidityFactor=1.8;
    if (W_UNIFORM()<0.5) {
      _BaselineAvailabilityToMosquitoes=0.2;
      _comorbidityFactor=0.2;
    }
  }  
  if (isOptionIncluded(modelVersion,tripleHet)) {
    _BaselineAvailabilityToMosquitoes=1.8;
    _comorbidityFactor=1.8;
    _treatmentSeekingFactor=0.2;
    if (W_UNIFORM()<0.5) {
      _BaselineAvailabilityToMosquitoes=0.2;
      _comorbidityFactor=0.2;
      _treatmentSeekingFactor=1.8;
    }
  }  
  if (isOptionIncluded(modelVersion, includesPKPD)) { 
    _drugs = list<Drug*>();
    _proxy = new DrugProxy(this);
  }
  _withinHostModel = new OldWithinHostModel(this);
}

Human::Human(istream& funit, CaseManagement* caseManagement, int simulationTime) 
  : _simulationTime(simulationTime), _caseManagement(caseManagement) {
  
  if (isOptionIncluded(modelVersion, includesPKPD)) { 
    _drugs = list<Drug*>();
    _proxy = new DrugProxy(this);
  }
  _withinHostModel = new OldWithinHostModel(this);
  // Reading human from file
  funit >> *this;
  _latestEvent.setCaseManagement(caseManagement);

  if ( _MOI <  0) {
    cout << "E MOI" << endl;
    exit(-3);
  }

}

Human::~Human(){
  if (_latestEvent.getTime() !=  missing_value){
    Simulation::gMainSummary->report(_latestEvent);
  }  
  clearAllInfections();
  if (isOptionIncluded(modelVersion, includesPKPD)) { // treatInfections() 
    list<Drug*>::const_iterator it;
    for(it=_drugs.begin(); it!=_drugs.end(); it++) {
      delete (*it);
    }
    _drugs.clear();
    delete _proxy;
  }
  delete _withinHostModel;
}

void Human::readFromFile(fstream& funit){

  funit >> *this;

  if ( _MOI <  0) {
    cout << "E MOI" << endl;
    exit(-3);
  }

}

void Human::updateInfection(double adultEIR, double availability, double expectedNumberOfInfections){

  introduceInfections(adultEIR, availability, expectedNumberOfInfections); 
  clearOldInfections();

  if ( mymodf(_simulationTime*interval, 5) ==  0) {
    for (int i=1;i<=3; i++) {
      _ylag[4-i]=_ylag[3-i];
    }
  }
  _ylag[0]=_totalDensity;
  _cumulativeYlag = _cumulativeY;

  _withinHostModel->calculateDensities();
}

void Human::treatInfections(){

  if (isOptionIncluded(modelVersion, includesPKPD)) { // treatInfections() 
    treatAllInfections();
    _proxy->decayDrugs();
  }

}

void Human::clearOldInfections(){

  int enddate;
  std::list<Infection*>::iterator iter=infections.begin();
  while(iter != infections.end()){
    enddate=(*iter)->getEndDate();
    //enddate=(*iter)->iData.startdate+(*iter)->iData.duration/interval;
    if (_simulationTime >= enddate) {
      delete *iter;
      iter=infections.erase(iter);
      _MOI=_MOI-1;
    }
    else{
      iter++;
    }
  }
}

void Human::clearAllInfections(){
  std::list<Infection*>::iterator i;
  for(i=infections.begin(); i != infections.end(); i++){
    //std::cout<<"clear all: "<<_ID<<std::endl;
    delete *i;
  }
  infections.clear();
  _MOI=0;
}

void Human::treatAllInfections(){
  std::list<Infection*>::iterator i;
  for(i=infections.begin(); i != infections.end(); i++){
    (*i)->setDensity((*i)->getDensity()*exp(-_proxy->calculateDrugsFactor((*i))));
  }
}

void Human::introduceInfections(double adultEIR, double availability, double expectedNumberOfInfections){
  double pInfectedstep;
  int Ninf;
  int i;
  //TODO: this code does not allow for variations in baseline availability
  //this is only likely to be relevant in some models but should not be
  //forgotten
  
  //Update pre-erythrocytic immunity
  if (isOptionIncluded(modelVersion,transHet) || isOptionIncluded(modelVersion,comorbTransHet) || isOptionIncluded(modelVersion,transTreatHet) || isOptionIncluded(modelVersion, tripleHet)) {
    _cumulativeEIRa+=(double)(interval*availability*adultEIR*(_BaselineAvailabilityToMosquitoes));
  }
  else {
    _cumulativeEIRa+=(double)interval*availability*adultEIR;
  }

  if ( expectedNumberOfInfections >  0.0000001) {
    Ninf=W_POISSON(expectedNumberOfInfections);
    if ( Ninf >  0) {
      for ( i=1;i<=Ninf; i++) {
        newInfection();
      }
    }
  }

  pInfectedstep=1-exp(-expectedNumberOfInfections);
  _pinfected=1.0-(1.0-pInfectedstep)*(1.0-_pinfected);
  _pinfected=std::min(_pinfected, 1.0);
  _pinfected=std::max(_pinfected, 0.0);

}

void Human::update(int simulationTime, TransmissionModel* transmissionModel){

  _simulationTime = simulationTime;
  double availabilityUpdate = transmissionModel->getRelativeAvailability(getAgeInYears());
  updateInterventionStatus(); 
  updateImmuneStatus();
  updateInfection(transmissionModel->calculateEIR(_simulationTime),
                        transmissionModel->getRelativeAvailability(getAgeInYears()),
                        transmissionModel->getExpectedNumberOfInfections(getPEVEfficacy(),
                                                                         getBaselineAvailabilityToMosquitoes(),
                                                                         getCumulativeEIRa(),
                                                                         availabilityUpdate*transmissionModel->calculateEIR(_simulationTime)));
  determineClinicalStatus();
  setWeight(transmissionModel->getWeight(getAgeInYears()));
  treatInfections();

}

void Human::setCaseManagement(CaseManagement* caseManagement){
  _caseManagement=caseManagement;
}

void Human::newInfection(){
   
  //std::cout<<"MOI "<<_MOI<<std::endl;
  if (_MOI <=  20) {
    _cumulativeInfections=_cumulativeInfections+1;
    infections.push_back(new Infection(_lastSPDose, _simulationTime));
    _MOI=_MOI+1;
  }
}

void Human::updateImmuneStatus(){
  if (immEffectorRemain < 1){
    _cumulativeh=(double)_cumulativeh*immEffectorRemain;
    _cumulativeY=(double)_cumulativeY*immEffectorRemain;
  }
  if (asexImmRemain < 1){
    _cumulativeh=(double)_cumulativeh*asexImmRemain/
      (1+(_cumulativeh*(1-asexImmRemain)/Infection::cumulativeHstar));
    _cumulativeY=(double)_cumulativeY*asexImmRemain/
      (1+(_cumulativeY*(1-asexImmRemain)/Infection::cumulativeYstar));
  }
}

void Human::doCM(int entrypoint){
 
  //TODO: implement age-specificity of drug dosing
  int decisionID;
  int nMedicate;
  int medicateID;
  int time;
  double qty;
  //char name;
  double ageyrs= getAgeInYears();
  decisionID=get_decision_id(entrypoint, ageyrs);
  nMedicate=get_n_medicate(decisionID, ageyrs);
  if ( nMedicate ==  0) {
    return ;
  }
  for ( medicateID=1;medicateID<=nMedicate; medicateID++) {
    qty=get_cmp_qty(decisionID, medicateID, ageyrs);
    time=get_cmp_time(decisionID, medicateID, ageyrs);
    //name=get_cmp_name(decisionID, medicateID, ageyrs);
    //call medicate(current,name,qty,time)
    cout << entrypoint;
    cout << ageyrs;
    //cout << name;
    cout << qty;
    cout << time << endl;
  }
}

void Human::medicate(string drugName, double qty, int time) {
  _proxy->medicate(drugName, qty, time);
}

void Human::determineClinicalStatus(){ //TODO: this function should not do case management
  //indirect death
  short iDeath;
  //indicates if this individual was treated successfully (ie parasites cleared)
  short drugEffect;
  // agegroup
  int agegrp;
  int fortnight;
  fortnight= round(14.0/interval);
  //TODOConversion
  //fortnight= (int)(14.0/interval);
  //fortnight=3;
  agegrp=ageGroup();
  //	Countdown to indirect mortality
  if ( _doomed <  0) {
    _doomed--;
  }
  iDeath=_latestEvent.indirectDeath(_simulationTime, _dateOfBirth, agegrp, _doomed);
  if (! iDeath) {
    drugEffect=this->defineEvent();
    if ( drugEffect) {
      if (! IPT) {
        if (! isOptionIncluded(modelVersion, includesPKPD)) {
          clearAllInfections();
        }
      }
      // TODO: maybe all this IPT stuff should be put in a different module
      if ( IPT) {
        if ( _latestEvent.getDiagnosis() ==  Diagnosis::SEVERE_MALARIA) {
          clearAllInfections();
        }
        else if(_simulationTime-_lastIptiOrPlacebo <=  fortnight) {
          clearAllInfections();
          // IPTi trials used quinine for fevers within 14 days of an ipti or placebo dose   
        }
        else if(_simulationTime-_lastSPDose <=  fortnight) {
          clearAllInfections();
          /*
            second line used if fever within 14 days of SP dose (ipti or treatment)
            TODO: if this code is to survive, then the iptiEffect values should be 
            symbolic constants
          */
        }
        else if( iptiEffect ==  2 ||  iptiEffect ==  12) {
          clearAllInfections();
          _lastSPDose=_simulationTime+1;
        }
        else if( iptiEffect ==  3 ||  iptiEffect ==  13) {
          clearAllInfections();
        }
        else if(iptiEffect >=  14 && iptiEffect < 30) {
          clearAllInfections();
        }
        else {
          _lastSPDose=_simulationTime+1;
          clearAllInfections();
          // SPAction will first act at the beginning of the next interval
        }
      }
    }
    // if drugEffect
  }
  // if not ideath
}

void Human::vaccinate(){
  
  double a;
  //Index to look up initial efficacy relevant for this dose.
  int effIndex;
  _lastVaccineDose=_lastVaccineDose+1;
  effIndex=std::min(_numberOfInitEff, _lastVaccineDose);
  if ( isOptionIncluded(vaccineType, preerythrocytic_reduces_h)) {
    if ( PEVInitialMeanEfficacy[effIndex - 1] <  1) {
      a=PEVefficacyB*PEVInitialMeanEfficacy[effIndex - 1]/(1.0-PEVInitialMeanEfficacy[effIndex - 1]);
      _PEVEfficacy=(W_BETA(a, (PEVefficacyB)));
    }
    else {
      _PEVEfficacy=1;
    }
  }
  if ( isOptionIncluded(vaccineType, erythrocytic_reduces_y)) {
    if ( BSVInitialMeanEfficacy[effIndex - 1] <  1) {
      a=BSVefficacyB*BSVInitialMeanEfficacy[effIndex - 1]/(1.0-BSVInitialMeanEfficacy[effIndex - 1]);
      _BSVEfficacy=(W_BETA(a, (BSVefficacyB)));
    }
    else {
      _BSVEfficacy=1;
    }
  }
  if ( isOptionIncluded(vaccineType, transmission_blocking_reduces_k)) {
    if ( TBVInitialMeanEfficacy[effIndex - 1] <  1) {
      a=TBVefficacyB*TBVInitialMeanEfficacy[effIndex - 1]/(1.0-TBVInitialMeanEfficacy[effIndex - 1]);
      _TBVEfficacy=(W_BETA(a, (TBVefficacyB)));
    }
    else {
      _TBVEfficacy=1;
    }
  }
}

void Human::updateInterventionStatus(){
   
  int agetstep;
  int nextDose;
  if ( vaccineType >  0) {
    agetstep=_simulationTime-_dateOfBirth;
    /*
      Update the effect of the vaccine
      We should assume the effect is maximal 25 days after vaccination
      TODO: consider the sensitivity of the predictions to the introduction 
      of a delay until the vaccine has reached max. efficacy.
    */
    if ( _lastVaccineDose >  0) {
      _PEVEfficacy=_PEVEfficacy*exp(-PEVdecay);
      _TBVEfficacy=_TBVEfficacy*exp(-TBVdecay);
      _BSVEfficacy=_BSVEfficacy*exp(-BSVdecay);
    }
    /*
      Determine eligibility for vaccination
      check for number of vaccine doses in the vaccinate subroutine
      TODO: The tstep conditional is appropriate if we assume there is no intervention during warmup
      It won't work if we introduce interventions into a scenario with a pre-existing intervention.
    */
    if ( Simulation::timeStep >  0) {
      nextDose=_lastVaccineDose+1;
      if (nextDose<=_numberOfEpiDoses){
        if ( (W_UNIFORM()) <  vaccineCoverage[nextDose - 1] &&  targetagetstep[nextDose - 1] ==  agetstep) {
          vaccinate();
          Simulation::gMainSummary->reportEPIVaccination(ageGroup());
        }
      }
    }
  }
  if ( ITN) {
    if ( (W_UNIFORM()) <  ITNdecay) {
      _ITN=false;/* bool */;
    }
  }
  if ( IPT) {
    if ( Simulation::timeStep >  0) {
      if (iptiEffect==14) {
        if ((modIntervalsPerYear(_simulationTime)-1)>=42 && (modIntervalsPerYear(_simulationTime)-1)<60) {
          // assumes 5-day intervals and Niakhar seasonality
          setLastSPDose();
        }
      }
      else if (iptiEffect==15) {
        if ((modIntervalsPerYear(_simulationTime)-1)>=48 && (modIntervalsPerYear(_simulationTime)-1)<66) {
          // assumes 5-day intervals and Niakhar seasonality
          setLastSPDose();
        }
      }
      else if (iptiEffect==16) {
        if ((modIntervalsPerYear(_simulationTime)-1)>=54 && (modIntervalsPerYear(_simulationTime)-1)<72) {
          // assumes 5-day intervals and Niakhar seasonality
          setLastSPDose();
        }
      }
      else if (iptiEffect==17) {
        if ((modIntervalsPerYear(_simulationTime)-1)>=60 && (modIntervalsPerYear(_simulationTime)-1)<78) {
          // assumes 5-day intervals and Niakhar seasonality
          setLastSPDose();
        }
      }
      else if (iptiEffect==18) {
        if ((modIntervalsPerYear(_simulationTime)-1)>=66 && (modIntervalsPerYear(_simulationTime)-1)<82) {
          // assumes 5-day intervals and Niakhar seasonality
          setLastSPDose();
        }
      }
     else if (iptiEffect==19) {
        if ((modIntervalsPerYear(_simulationTime)-1)>=36 && (modIntervalsPerYear(_simulationTime)-1)<54) {
          // assumes 5-day intervals and Niakhar seasonality
          setLastSPDose();
        }
      }
     else if (iptiEffect==20) {
        if ((modIntervalsPerYear(_simulationTime)-1)>=30 && (modIntervalsPerYear(_simulationTime)-1)<48) {
          // assumes 5-day intervals and Niakhar seasonality
          setLastSPDose();
        }
      }
     else if (iptiEffect==21) {
        if ((modIntervalsPerYear(_simulationTime)-1)>=24 && (modIntervalsPerYear(_simulationTime)-1)<42) {
          // assumes 5-day intervals and Niakhar seasonality
          setLastSPDose();
        }
      }
     else if (iptiEffect==22) {
        if ((modIntervalsPerYear(_simulationTime)-1)>=18 && (modIntervalsPerYear(_simulationTime)-1)<42) {
          // assumes 5-day intervals and Niakhar seasonality
          setLastSPDose();
        }
      }
      else {
        setLastSPDose();
      }
    }
  }
}

void Human::setLastSPDose() {
  double rnum;
  int agetstep;
  int i;
  agetstep=_simulationTime-_dateOfBirth;
  for ( i=1;i<=numberOfIPTiDoses; i++) {
    if ( iptiTargetagetstep[i - 1] == agetstep) {
      rnum=W_UNIFORM();
      if ( (rnum) <  iptiCoverage[i - 1]) {
        _lastIptiOrPlacebo=_simulationTime;
        /*
          iptiEffect denotes treatment or placebo group
          and also the treatment given when sick (trial-dependent)
        */
        if (iptiEffect >=  10) {
          _lastSPDose=_simulationTime;
          Simulation::gMainSummary->reportIPTDose(ageGroup());					 
        }
      }
    }
  }
}

void Human::summarize(){
    
  double age = getAgeInYears();
  if (get_is_ipti() == 1 && _simulationTime-_tLastTreatment >= 1 && _simulationTime-_tLastTreatment <=  4) {
    return ;
  }
  Simulation::gMainSummary->addToHost(age,1);
  if ( _MOI >  0) {
    Simulation::gMainSummary->addToInfectedHost(age,1);
    Simulation::gMainSummary->addToTotalInfections(age, _MOI);
    Simulation::gMainSummary->addToTotalPatentInfections(age, _patentinfections);

  }
  if ( _totalDensity >  detectionlimit) {
    Simulation::gMainSummary->addToPatentHost(age, 1);
    Simulation::gMainSummary->addToSumLogDensity(age, log(_totalDensity));
  }
  Simulation::gMainSummary->addToExpectedInfected(age, _pinfected);
  Simulation::gMainSummary->addToPyrogenicThreshold(age, _pyrogenThres);
  Simulation::gMainSummary->addToSumX(age, log(_pyrogenThres+1.0));
}

void initDrugAction () {
  //Part of this will be auto generated

}


int Human::ageGroup() const{
  double ageyears=getAgeInYears();
  return Simulation::gMainSummary->ageGroup(ageyears);
}

double Human::getAgeInYears() const{
  return 1.0*((_simulationTime-_dateOfBirth)*interval)/daysInYear;
}

double Human::getAgeInYears(int time) const{
  return 1.0*((time-_dateOfBirth)*interval)/daysInYear;
}

double Human::Ystar(){
  int i;
  //Number of categories in the numerical approx. below
  const int n= 11;
  const double delt= 1.0/n;
  double valYstar=_pyrogenThres;
  //Numerical approximation to equation 2, AJTMH p.57
  for ( i=1;i<=n; i++) {
    valYstar=valYstar+_totalDensity*alpha14*interval*delt/((Ystar1_26+_totalDensity)*(Ystar2_13+valYstar))-smuY*valYstar*delt;
  }
  return valYstar;
}

double Human::infectiousness(){
  double transmit;
  double x;
  int agetstep;
  //Infectiousness parameters: see AJTMH p.33, tau=1/sigmag**2 
  static const double beta2=0.46;
  static const double beta3=0.17;
  static const double tau= 0.066;
  static const double mu= -8.1;
  double pone;
  double zval;
  double valinfectiousness;
  agetstep=_simulationTime-_dateOfBirth;
  /*
    Original infectiousness model based on 5 day intervals updates
    lagged variables only every 5 days and cannot compute infectiousness
    for the first 20 days of the simulation
  */
  if (( agetstep*interval >  20) && ( _simulationTime*interval >  20)) {
    x=_ylag[1]+beta2*_ylag[2]+beta3*_ylag[3];
    if ( x <  0.001) {
      transmit=0.0;
    }
    else {
      zval=(log(x)+mu)/sqrt(1.0/tau);
      pone=W_UGAUSS_P((zval));
      transmit=(pone*pone);
    }
  }
  else {
    transmit=0.0;
  }
  //transmit has to be between 0 and 1
  transmit=std::max(transmit, 0.0);
  transmit=std::min(transmit, 1.0);
  //	Include here the effect of transmission-blocking vaccination
  valinfectiousness=transmit*(1-_TBVEfficacy);
  return valinfectiousness;
}

bool Human::uncomplicatedEvent(bool isMalaria){
  int nextRegimen;
  int entrypoint;
  int agegroup;
  bool returnValue=false;
  //ageGroup is not optimized
  agegroup=ageGroup();
  if ( isMalaria) {
    entrypoint=Diagnosis::UNCOMPLICATED_MALARIA;
  }
  else {
    entrypoint=Diagnosis::NON_MALARIA_FEVER;
  }
  if ( isOptionIncluded(modelVersion, caseManagementV2)) {
    returnValue=true;
  }
  else {
    nextRegimen=_caseManagement->getNextRegimen(_simulationTime, entrypoint, _tLastTreatment, _latestRegimen);
    //call doCM(nextRegimen,0.1)
    if (_caseManagement->getProbabilityGetsTreatment(nextRegimen-1)*_treatmentSeekingFactor > (W_UNIFORM())){
      if (isOptionIncluded(modelVersion, includesPKPD)){
        _latestEvent.update(_simulationTime, agegroup, entrypoint, Outcome::PARASITES_PKPD_DEPENDENT_RECOVERS_OUTPATIENTS);
        /*
          TODO: uncomplicatedEvent forces a call of pk PD model
          in the event that there is no treatment it should remain .false.
          TODO: 
          call medicate(currInd, "cq ", 300.0, 0)
        */
        returnValue=true;
      }
      else {
        //PostTiagoBUGFIX: remove order of or in order to prevent short circuit
        //if (isnan(pParasitesCleared[nextRegimen-1]) ||( pParasitesCleared[nextRegimen - 1] >  (W_UNIFORM()))) {
        if ((_caseManagement->getProbabilityParasitesCleared(nextRegimen-1) > W_UNIFORM())  || 
            (isnan(_caseManagement->getProbabilityParasitesCleared(nextRegimen-1)))){
          _latestEvent.update(_simulationTime, agegroup, entrypoint, Outcome::PARASITES_ARE_CLEARED_PATIENT_RECOVERS_OUTPATIENTS);
          returnValue=true;         
        }
        else {
          _latestEvent.update(_simulationTime, agegroup, entrypoint, Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_OUTPATIENTS);
        }
      }
      _latestRegimen=nextRegimen;
      _tLastTreatment=_simulationTime;
      Simulation::gMainSummary->reportTreatment(agegroup, _latestRegimen);
    }
    else {
      _latestEvent.update(_simulationTime, agegroup, entrypoint, Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_NON_TREATED);
    }
  }
  return returnValue;
}

bool Human::severeMalaria(){
  /*
    DOCU
    Set doomed=4 if the patient dies.
  */

  int nextRegimen;
  double q[9];
  double p2;
  double p3;
  double p4;
  double p5;
  double p6;
  double p7;
  int isAdult;
  double prandom;
  bool returnValue=false;
  int agegroup=ageGroup();
  prandom=(W_UNIFORM());
  isAdult=2;
  if (getAgeInYears() >= 5.0) {
    isAdult=1;
  }
  if ( isOptionIncluded(modelVersion, caseManagementV2)) {
    returnValue=true;      
  }
  /*
    TODO: is there a reason we cannot first decide if the patient is treated, then increase latestRegimen
    and latestTreatment, instead of resetting it if not treated? (Can't think of one, but
    do we want to change this section of code rather than just introducing the new alternative (TS))
  */ 
  nextRegimen=_caseManagement->getNextRegimen(_simulationTime, Diagnosis::SEVERE_MALARIA, _tLastTreatment, _latestRegimen);
  //call doCM(nextRegimen,0.1)
  p2=_caseManagement->getProbabilityGetsTreatment(nextRegimen-1)*_treatmentSeekingFactor;
  p3=_caseManagement->getCureRate(nextRegimen-1);
  p4=_caseManagement->caseFatality(getAgeInYears());
  /*
    p4 is the hospital case-fatality rate from Tanzania
    p5 here is the community threshold case-fatality rate
  */
  p5=_caseManagement->getCommunityCaseFatalityRate(p4);
  p6=_caseManagement->getProbabilitySequelaeTreated(isAdult-1);
  p7=_caseManagement->getProbabilitySequelaeUntreated(isAdult-1);
  //	Community deaths
  q[0]=(1-p2)*p5;
  //	Community sequelae
  q[1]=q[0]+(1-p2)*(1-p5)*p7;
  //	Community survival
  q[2]=q[1]+(1-p2)*(1-p5)*(1-p7);
  //	Parasitological failure deaths
  q[3]=q[2]+p2*p5*(1-p3);
  //	Parasitological failure sequelae
  q[4]=q[3]+p2*(1-p3)*(1-p5)*p7;
  //	Parasitological failure survivors
  q[5]=q[4]+p2*(1-p3)*(1-p5)*(1-p7);
  //	Parasitological success deaths
  q[6]=q[5]+p2*p3*p4;
  //	Parasitological success sequelae
  q[7]=q[6]+p2*p3*(1-p4)*p6;
  //	Parasitological success survival
  q[8]=q[7]+p2*p3*(1-p4)*(1-p6);
  /*
    if (q(5).lt.1) stop
    NOT TREATED
  */
  if ( q[0] >  prandom) {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PATIENT_DIES_NON_TREATED);
    _doomed  = 4;
  }
  else if( q[1] >  prandom) {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_NOT_CLEARED_PATIENT_HAS_SEQUELAE_NON_TREATED);
  }
  else if( q[2] >  prandom) {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_NON_TREATED);
  }
  else if( q[3] >  prandom) {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PATIENT_DIES_INPATIENTS);
    _tLastTreatment = _simulationTime;
    _latestRegimen = nextRegimen;
    Simulation::gMainSummary->reportTreatment(agegroup,_latestRegimen);
    _doomed  = 4;
  }
  else if( q[4] >  prandom) {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_NOT_CLEARED_PATIENT_HAS_SEQUELAE_INPATIENTS);
    _tLastTreatment = _simulationTime;
    _latestRegimen = nextRegimen;
    Simulation::gMainSummary->reportTreatment(agegroup,_latestRegimen);
  }
  else if( q[5] >  prandom) {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_INPATIENTS);
    _tLastTreatment = _simulationTime;
    _latestRegimen = nextRegimen;
    Simulation::gMainSummary->reportTreatment(agegroup,_latestRegimen);
  }
  else if( q[6] >  prandom) {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PATIENT_DIES_INPATIENTS);
    _tLastTreatment = _simulationTime;
    _latestRegimen = nextRegimen;
    returnValue=true;
    _doomed  = 4;
    Simulation::gMainSummary->reportTreatment(agegroup,_latestRegimen);
  }
  else if( q[7] >  prandom) {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_ARE_CLEARED_PATIENT_HAS_SEQUELAE_INPATIENTS);
    _tLastTreatment = _simulationTime;
    _latestRegimen = nextRegimen;
    returnValue=true;
    Simulation::gMainSummary->reportTreatment(agegroup,_latestRegimen);
  }
  else if( q[8] >=  prandom) {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_ARE_CLEARED_PATIENT_RECOVERS_INPATIENTS);
    _tLastTreatment = _simulationTime;
    _latestRegimen = nextRegimen;
    returnValue=true;      
    Simulation::gMainSummary->reportTreatment(agegroup,_latestRegimen);
  }
  else {
    cout << "SM: ONE" << endl;

    cout << prandom;
    cout << q[0];
    cout << q[1];
    cout << q[2];
    cout << q[3];
    cout << q[4];
    cout << q[5];
    cout << q[6];
    cout << q[7];
    cout << q[8] << endl;
    
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_ARE_CLEARED_PATIENT_RECOVERS_INPATIENTS);
    _tLastTreatment = _simulationTime;
    _latestRegimen = nextRegimen;
    returnValue=true;        
    Simulation::gMainSummary->reportTreatment(agegroup,_latestRegimen);
  }
  return returnValue;
}

bool Human::defineEvent(){
  //double q0;
  double prSevereEpisode;
  double prEpisode;
  double IncidenceDensity;
  double prNonMalariaFever;
  double indirectRisk;
  double pCoinfection;
  double severeMalThreshold;
  double ageYears;
  bool valdefineEvent =false;
  ageYears=getAgeInYears();
  if ( isOptionIncluded(modelVersion, predeterminedEpisodes)) {
    prEpisode=0;
    if ( _timeStepMaxDensity > _pyrogenThres) {
      prEpisode=1;
    }
  }
  else {
    if ( isOptionIncluded(modelVersion, MuellerMorbidityModel)) {
      IncidenceDensity=rateMultiplier_31*(pow(_totalDensity, densityExponent_32))/(1.0*intervalsPerYear);
      prEpisode=1-exp(-IncidenceDensity);
    }
    else {
      prEpisode=1-1/(1+(_timeStepMaxDensity/_pyrogenThres));
    }
  }
  //Fixed severe threshold
  severeMalThreshold=sevMal_21+1;
  prSevereEpisode=1-1/(1+_timeStepMaxDensity/severeMalThreshold);
  //Decide whether a clinical episode occurs and if so, which type
  pCoinfection=comorbintercept_24/(1+ageYears/critAgeComorb_30);
  pCoinfection=pCoinfection*_comorbidityFactor;
  if ((W_UNIFORM())< prEpisode) {
    if ((W_UNIFORM())< prSevereEpisode){
      valdefineEvent=severeMalaria();
    }
    else if( (W_UNIFORM()) < pCoinfection) {
      valdefineEvent=severeMalaria();
    }
    else {
      valdefineEvent=uncomplicatedEvent(true);
    }
    /*
      Indirect mortality	
      IndirectRisk is the probability of dying from indirect effects of malaria
      conditional on not having an acute attack of malaria
    */
    indirectRisk=indirRiskCoFactor_18/(1+ageYears/critAgeComorb_30);
    indirectRisk=indirectRisk*_comorbidityFactor;
    if ( (W_UNIFORM()) < indirectRisk) {
      if (_doomed ==  0) {
        _doomed=-1;
      }
    }
    // Penalizing immunity for clinical episodes	
    if ( isOptionIncluded(modelVersion, penalisationEpisodes)) {
      _cumulativeY=(double)_cumulativeYlag+(-(immPenalty_22*(_cumulativeY-_cumulativeYlag)));
      if (_cumulativeY <  0) {
        _cumulativeY=0.0;
      }
    }
  }
  else if( isOptionIncluded(modelVersion, nonMalariaFevers)) {
    prNonMalariaFever=pCoinfection*RelativeRiskNonMalariaFever;
    if ((W_UNIFORM()) < prNonMalariaFever){
      valdefineEvent=uncomplicatedEvent(false);
    }
  }
  return valdefineEvent;
}


double Human::calculateSelectionCoefficient (Infection inf) {
  return 1.0;
}

list<Drug*>* Human::getDrugs() {
  return &_drugs;
}

ostream& operator<<(ostream& out, const Human& human){

  out << human._cumulativeInfections << endl; 
  out << human._dateOfBirth << endl; 
  out << human._doomed << endl; 
  out << human._ID << endl ; 
  out << human._lastVaccineDose << endl;
  out << human._latestRegimen << endl; 
  out << human._tLastTreatment << endl; 
  out << human._MOI << endl; 
  out << human._patentinfections << endl; 
  out << human._BSVEfficacy << endl; 
  out << human._cumulativeYlag << endl; 
  out << human._cumulativeEIRa << endl; 
  out << human._cumulativeh << endl; 
  out << human._cumulativeY << endl; 
  out << human._timeStepMaxDensity << endl; 
  out << human._innateImmunity << endl; 
  out << human._PEVEfficacy << endl; 
  out << human._pinfected << endl; 
  out << human._ptransmit << endl; 
  out << human._pyrogenThres << endl; 
  out << human._TBVEfficacy << endl; 
  out << human._totalDensity << endl; 
  out << human._ylag[0] << endl; 
  out << human._ylag[1] << endl; 
  out << human._ylag[2] << endl; 
  out << human._ylag[3] << endl; 
  out << boolalpha << human._ITN << endl; 
  out << human._lastSPDose << endl; 
  out << human._lastIptiOrPlacebo << endl; 
  out << human._comorbidityFactor << endl;  
  out << human._treatmentSeekingFactor << endl; 
  out << human._BaselineAvailabilityToMosquitoes << endl; 
  out << human._latestEvent;   
  out << *human._withinHostModel;

    for(std::list<Infection*>::const_iterator iter=human.infections.begin(); iter != human.infections.end(); iter++)
      out << **iter;
     
  if (isOptionIncluded(modelVersion, includesPKPD)) {
    out << *(human._proxy);
    out << human._drugs.size() << endl;
    list<Drug*>::const_iterator it;
    for (it=human._drugs.begin(); it!=human._drugs.end(); it++) {
      out << **it;
    }
  }
  return out;

}

istream& operator>>(istream& in, Human& human){

  in >> human._cumulativeInfections; 
  in >> human._dateOfBirth; 
  in >> human._doomed; 
  in >> human._ID; 
  in >> human._lastVaccineDose; 
  in >> human._latestRegimen; 
  in >> human._tLastTreatment; 
  in >> human._MOI; 
  in >> human._patentinfections; 
  in >> human._BSVEfficacy; 
  in >> human._cumulativeYlag; 
  in >> human._cumulativeEIRa; 
  in >> human._cumulativeh; 
  in >> human._cumulativeY; 
  in >> human._timeStepMaxDensity; 
  in >> human._innateImmunity; 
  in >> human._PEVEfficacy; 
  in >> human._pinfected; 
  in >> human._ptransmit; 
  in >> human._pyrogenThres; 
  in >> human._TBVEfficacy; 
  in >> human._totalDensity; 
  in >> human._ylag[0]; 
  in >> human._ylag[1]; 
  in >> human._ylag[2]; 
  in >> human._ylag[3]; 
  in >> boolalpha >> human._ITN; 
  in >> human._lastSPDose; 
  in >> human._lastIptiOrPlacebo; 
  in >> human._comorbidityFactor;  
  in >> human._treatmentSeekingFactor; 
  in >> human._BaselineAvailabilityToMosquitoes; 
  in >> human._latestEvent;
  in >> *(human._withinHostModel);

  for(int i=0;i<human._MOI;++i) {
    Infection* infection = new Infection();
    in >> *infection;
    human.infections.push_back(infection);
  }

  if (isOptionIncluded(modelVersion, includesPKPD)) {
    int numDrugs;
    in >> *(human._proxy);
    in >> numDrugs;
    for (int i=0; i<numDrugs; i++) {
      Drug* drug = new Drug("", "", 0, 0);
      in >> *drug;
      human._drugs.push_back(drug);
    }
    
  }

  return in;

}


