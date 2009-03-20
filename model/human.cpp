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
#include "DescriptiveInfection.h"


/*
  Constants common to all humans
  Decay in anti-toxic immunity
*/
//TODOConversion
#ifdef _WIN32
#define isnan(x) ((x) != (x))
#endif


double Human::sigma_i;
double Human::sevMal_21;
double Human::critAgeComorb_30;
double Human::comorbintercept_24;
double Human::indirRiskCoFactor_18;
double Human::immPenalty_22;
double Human::asexImmRemain;
double Human::immEffectorRemain;
double Human::BaselineAvailabilityShapeParam;
double Human::detectionlimit;
double Human::baseEntoAvailability;
double Human::baseProbMosqSurvivalBiting;
double Human::baseProbMosqSurvivalResting;

const double Human::wtprop[nwtgrps] = { 0.116547265, 0.152531009, 0.181214575, 0.202146126, 0.217216287, 0.237405732, 0.257016899, 0.279053187, 0.293361286, 0.309949502, 0.334474135, 0.350044993, 0.371144279, 0.389814144, 0.412366341, 0.453, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };

void Human::initHumanParameters () {	// static
  /*
    Init parameters that are common to all humans
  */
  
  baseEntoAvailability = 1.0;		// FIXME: get from xml
  baseProbMosqSurvivalBiting = 1.0;	// FIXME: get from xml
  baseProbMosqSurvivalResting = 1.0;	// FIXME: get from xml
  double densitybias;
  sigma_i=sqrt(getParameter(Params::SIGMA_I_SQ));
  indirRiskCoFactor_18=(1-exp(-getParameter(Params::INDIRECT_RISK_COFACTOR)));
  sevMal_21=getParameter(Params::SEVERE_MALARIA_THRESHHOLD);
  immPenalty_22=1-exp(getParameter(Params::IMMUNITY_PENALTY));
  comorbintercept_24=1-exp(-getParameter(Params::COMORBIDITY_INTERCEPT));
  immEffectorRemain=exp(-getParameter(Params::IMMUNE_EFFECTOR_DECAY));
  asexImmRemain=exp(-getParameter(Params::ASEXUAL_IMMUNITY_DECAY));
  critAgeComorb_30=getParameter(Params::CRITICAL_AGE_FOR_COMORBIDITY);
  /*
    TODO: This densitiybias function should be part of the scenario description XML, not the parameter element.
    or maybe it should be a parameter, as we want to fit it... but the garki analysis numbers are a bit dangerous
    add an attribute to scenario.xml densityQuantification="malariaTherapy|garki|other"
  */
  if (( get_analysis_no() <  22) || ( get_analysis_no() >  30)) {
    densitybias=getParameter(Params::DENSITY_BIAS_NON_GARKI);
  }
  else {
    densitybias=getParameter(Params::DENSITY_BIAS_GARKI);
  }
  detectionlimit=get_detectionlimit()*densitybias;
  BaselineAvailabilityShapeParam=getParameter(Params::BASELINE_AVAILABILITY_SHAPE);
}

// Testing Ctor
Human::Human() {
  _withinHostModel = new OldWithinHostModel(this);
  _morbidityModel=MorbidityModel::createMorbidityModel();
  if (Global::modelVersion & INCLUDES_PK_PD) {
    _drugs = list<Drug*>();
    _proxy = new DrugProxy(this);
  }
}

// Create new human
Human::Human(int ID, int dateOfBirth, CaseManagementModel* caseManagement, int simulationTime) 
  : _simulationTime(simulationTime), _caseManagement(caseManagement)
{
  //FIXME: should be partially random:
  _entoAvailability = baseEntoAvailability;
  _probMosqSurvivalBiting = baseProbMosqSurvivalBiting;
  _probMosqSurvivalResting = baseProbMosqSurvivalResting;
  
  //std::cout<<"newH:ID dateOfBirth "<<ID<<" "<<dateOfBirth<<std::endl;
  _BSVEfficacy=0.0;
  _cumulativeEIRa=0.0;
  _cumulativeh=0.0;
  _cumulativeInfections=0;
  _cumulativeY=0.0;
  _cumulativeYlag=0.0;
  _dateOfBirth=dateOfBirth;
  if (_dateOfBirth > simulationTime) {
    // This test may be totally unnecessary; it was done in oldWithinHostModel
    cerr << "date of birth in future!" << endl;
    throw 0;
  }
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
  _TBVEfficacy=0.0;
  _totalDensity=0.0;
  for (int i=1;i<=4; i++) {
    _ylag[i-1]=0.0;
  }
  _latestEvent.setTime(missing_value);
  _latestEvent.setCaseManagement(caseManagement);
  _lastSPDose=missing_value;
  _lastIptiOrPlacebo=missing_value;
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    _BaselineAvailabilityToMosquitoes=(W_GAMMA((BaselineAvailabilityShapeParam), (BaselineAvailabilityMean/BaselineAvailabilityShapeParam)));
  }
  else if(Global::modelVersion & LOGNORMAL_MASS_ACTION) {
    _BaselineAvailabilityToMosquitoes=(W_LOGNORMAL((log(BaselineAvailabilityMean))-(0.5*pow(BaselineAvailabilityShapeParam, 2)), (BaselineAvailabilityShapeParam)));
  }
  else if (Global::modelVersion & TRANS_HET) {
    _BaselineAvailabilityToMosquitoes=0.2;
    if (W_UNIFORM() < 0.5) {            
      _BaselineAvailabilityToMosquitoes=1.8;
    }
  }
  else {
    _BaselineAvailabilityToMosquitoes=BaselineAvailabilityMean;
  }
  if (Global::modelVersion & COMORB_HET) {
    _comorbidityFactor=0.2;
    if (W_UNIFORM() < 0.5) {            
      _comorbidityFactor=1.8;
    }	
  }
  else {
    _comorbidityFactor=1.0;
  }
  if (Global::modelVersion & TREAT_HET) {
    _treatmentSeekingFactor=0.2;
    if (W_UNIFORM() < 0.5) {            
      _treatmentSeekingFactor=1.8;
    }	
  }
  else {
    _treatmentSeekingFactor=1.0;
  }
  if (Global::modelVersion & TRANS_TREAT_HET) {
    _treatmentSeekingFactor=0.2;
    _BaselineAvailabilityToMosquitoes=1.8;
    if (W_UNIFORM()<0.5) {
      _treatmentSeekingFactor=1.8;
      _BaselineAvailabilityToMosquitoes=0.2;
    }
  }
  if (Global::modelVersion & COMORB_TREAT_HET) {
    _comorbidityFactor=0.2;
    _treatmentSeekingFactor=1.8;
    if (W_UNIFORM()<0.5) {
      _treatmentSeekingFactor=0.2;
      _comorbidityFactor=1.8;
    }
  }
  if (Global::modelVersion & COMORB_TRANS_HET) {
    _BaselineAvailabilityToMosquitoes=1.8;
    _comorbidityFactor=1.8;
    if (W_UNIFORM()<0.5) {
      _BaselineAvailabilityToMosquitoes=0.2;
      _comorbidityFactor=0.2;
    }
  }  
  if (Global::modelVersion & TRIPLE_HET) {
    _BaselineAvailabilityToMosquitoes=1.8;
    _comorbidityFactor=1.8;
    _treatmentSeekingFactor=0.2;
    if (W_UNIFORM()<0.5) {
      _BaselineAvailabilityToMosquitoes=0.2;
      _comorbidityFactor=0.2;
      _treatmentSeekingFactor=1.8;
    }
  }  
  if (Global::modelVersion & INCLUDES_PK_PD) { 
    _drugs = list<Drug*>();
    _proxy = new DrugProxy(this);
  }
  _withinHostModel = new OldWithinHostModel(this);
  _morbidityModel=MorbidityModel::createMorbidityModel();
}

// Load human from checkpoint
Human::Human(istream& funit, CaseManagementModel* caseManagement, int simulationTime) 
  : _simulationTime(simulationTime), _caseManagement(caseManagement)
{
  if (Global::modelVersion & INCLUDES_PK_PD) { 
    _drugs = list<Drug*>();
    _proxy = new DrugProxy(this);
  }
  _withinHostModel = new OldWithinHostModel(this);
  _morbidityModel=MorbidityModel::createMorbidityModel();
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
  if (Global::modelVersion & INCLUDES_PK_PD) { // treatInfections() 
    list<Drug*>::const_iterator it;
    for(it=_drugs.begin(); it!=_drugs.end(); it++) {
      delete (*it);
    }
    _drugs.clear();
    delete _proxy;
  }
  delete _withinHostModel;
  delete _morbidityModel;
}

void Human::readFromFile(fstream& funit){

  funit >> *this;

  if ( _MOI <  0) {
    cout << "E MOI" << endl;
    exit(-3);
  }

}

void Human::updateInfection(double expectedInfectionRate, double expectedNumberOfInfections){

  introduceInfections(expectedInfectionRate, expectedNumberOfInfections); 
  clearOldInfections();

  if ((_simulationTime*Global::interval) % 5 ==  0) {
    for (int i=1;i<=3; i++) {
      _ylag[4-i]=_ylag[3-i];
    }
  }
  _ylag[0]=_totalDensity;
  _cumulativeYlag = _cumulativeY;

  _withinHostModel->calculateDensities();
}

void Human::treatInfections(){

  if (Global::modelVersion & INCLUDES_PK_PD) { // treatInfections() 
    treatAllInfections();
    _proxy->decayDrugs();
  }

}

void Human::clearOldInfections(){

  int enddate;
  std::list<Infection*>::iterator iter=infections.begin();
  while(iter != infections.end()){
    enddate=(*iter)->getEndDate();
    //enddate=(*iter)->iData.startdate+(*iter)->iData.duration/Global::interval;
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

void Human::introduceInfections(double expectedInfectionRate, double expectedNumberOfInfections){
  //TODO: this code does not allow for variations in baseline availability
  //this is only likely to be relevant in some models but should not be
  //forgotten
  
  //Update pre-erythrocytic immunity
  if (Global::modelVersion & 
      (TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET)) {
    _cumulativeEIRa+=(double)(Global::interval*expectedInfectionRate*(_BaselineAvailabilityToMosquitoes));
  }
  else {
    _cumulativeEIRa+=(double)Global::interval*expectedInfectionRate;
  }

  if (expectedNumberOfInfections > 0.0000001) {
    int Ninf = W_POISSON(expectedNumberOfInfections);
    if (Ninf > 0) {
      for (int i=1;i<=Ninf; i++) {
        newInfection();
      }
    }
  }

  double pInfectedstep = 1.0 - exp(-expectedNumberOfInfections);
  _pinfected = 1.0 - (1.0-pInfectedstep) * (1.0-_pinfected);
  _pinfected = std::min(_pinfected, 1.0);
  _pinfected = std::max(_pinfected, 0.0);
}

void Human::update(int simulationTime, TransmissionModel* transmissionModel) {
  _simulationTime = simulationTime;
  double expectedInfectionRate = transmissionModel->getRelativeAvailability(getAgeInYears()) * transmissionModel->calculateEIR(_simulationTime, *this);
  
  updateInterventionStatus(); 
  updateImmuneStatus();
  updateInfection(expectedInfectionRate,
                  transmissionModel->getExpectedNumberOfInfections(*this, expectedInfectionRate));
  determineClinicalStatus();
  _weight = 120.0 * wtprop[TransmissionModel::getAgeGroup(getAgeInYears())];
  treatInfections();
}

void Human::setCaseManagement(CaseManagementModel* caseManagement) {
  _caseManagement=caseManagement;
}

void Human::newInfection(){
  //std::cout<<"MOI "<<_MOI<<std::endl;
  if (_MOI <=  20) {
    _cumulativeInfections=_cumulativeInfections+1;
    infections.push_back(new DescriptiveInfection(_lastSPDose, _simulationTime));
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
  //TODO: This is a rough and quick implementation, which could perhaps be improved.
  
  double ageyrs = getAgeInYears();
  if (getCaseManagements() == NULL) return;
  const CaseManagements::CaseManagementSequence managements = getCaseManagements()->getCaseManagement();
  if (managements.size() == 0) return;
  const CaseManagement* caseManagement = NULL;
  for (CaseManagements::CaseManagementConstIterator it = managements.begin(); it != managements.end(); ++it)
    if (ageyrs < it->getMaxAgeYrs() &&
        (!it->getMinAgeYrs().present() || it->getMinAgeYrs().get() <= ageyrs))
      caseManagement = &*it;
  if (caseManagement == NULL) {
    cerr << "No case management for age " << ageyrs << endl;
    throw 0;
  }
  
  const CaseType::EndPointSequence* caseTypeSeq;
  if (entrypoint == 1) {
    caseTypeSeq = &caseManagement->getUc1().getEndPoint();
  } else if (entrypoint == 2) {
    caseTypeSeq = &caseManagement->getUc2().getEndPoint();
  } else if (entrypoint == 3) {
    caseTypeSeq = &caseManagement->getSev().getEndPoint();
  } else if (entrypoint == 4) {
    caseTypeSeq = &caseManagement->getNmf().getEndPoint();
  } else {
    cerr << "invalid entrypoint" << endl;
    throw 0;
  }
  double randCum = W_UNIFORM();
  int decisionID = -1;
  for (CaseType::EndPointConstIterator it = caseTypeSeq->begin(); it != caseTypeSeq->end(); ++it) {
    randCum -= it->getP();
    if (randCum < 0) {
      decisionID = it->getDecision();
      break;
    }
  }
  if (decisionID < 0) {
    cerr << "Sum of probabilities of case management end-points for some severity type less than 1" << endl;
    throw 0;
  }
  
  //FIXME: build list of decisions by ID and use
  const Decisions::DecisionSequence& decisions = caseManagement->getDecisions().getDecision();
  if ((int)decisions.size() < decisionID) {
    cerr << "A decision for a case-management end-point doesn't exist (number "<<decisionID<<")!" << endl;
    return;
  }
  const Decision::MedicateSequence& medicates = decisions[decisionID-1].getMedicate();
  if (medicates.size() == 0)
    return;
  
  for (size_t medicateID=0;medicateID<medicates.size(); medicateID++) {
    double qty=medicates[medicateID].getQty();
    int time=medicates[medicateID].getTime();
    string name=medicates[medicateID].getName();
    medicate(name,qty,time);
  }
}

void Human::medicate(string drugName, double qty, int time) {
  _proxy->medicate(drugName, qty, time);
}

void Human::determineClinicalStatus(){ //TODO: this function should not do case management
  int fortnight = int((14.0/Global::interval)+0.5);	// round to nearest
  //TODOConversion
  //fortnight= (int)(14.0/Global::interval);
  //fortnight=3;
  //	Countdown to indirect mortality
  if ( _doomed <  0) {
    _doomed--;
  }
  //indirect death
  bool iDeath=_latestEvent.indirectDeath(_simulationTime, _dateOfBirth, ageGroup(), _doomed);
  if (! iDeath) {
  //indicates if this individual was treated successfully (ie parasites cleared)
    bool drugEffect=this->defineEvent();
    if ( drugEffect) {
      if (!IPTIntervention::IPT) {
        if (!(Global::modelVersion & INCLUDES_PK_PD)) {
          clearAllInfections();
        }
      }
      // TODO: maybe all this IPT stuff should be put in a different module
      if (IPTIntervention::IPT) {
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
# define iptiEffect IPTIntervention::iptiEffect
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
          // SPAction will first act at the beginning of the next Global::interval
        }
      }
    }
    // if drugEffect
  }
  // if not ideath
}

void Human::vaccinate(){
  //Index to look up initial efficacy relevant for this dose.
  if (Vaccine::PEV.active)
    _PEVEfficacy = Vaccine::PEV.getEfficacy(_lastVaccineDose);
  
  if (Vaccine::BSV.active)
    _BSVEfficacy = Vaccine::BSV.getEfficacy(_lastVaccineDose);
  
  if (Vaccine::TBV.active)
    _TBVEfficacy = Vaccine::TBV.getEfficacy(_lastVaccineDose);
  
  ++_lastVaccineDose;
}

void Human::updateInterventionStatus() {
  int agetstep;
  int nextDose;
  if (Vaccine::anyVaccine) {
    agetstep=_simulationTime-_dateOfBirth;
    /*
      Update the effect of the vaccine
      We should assume the effect is maximal 25 days after vaccination
      TODO: consider the sensitivity of the predictions to the introduction 
      of a delay until the vaccine has reached max. efficacy.
    */
    if ( _lastVaccineDose >  0) {
      _PEVEfficacy *= Vaccine::PEV.decay;
      _TBVEfficacy *= Vaccine::TBV.decay;
      _BSVEfficacy *= Vaccine::BSV.decay;
    }
    /*
      Determine eligibility for vaccination
      check for number of vaccine doses in the vaccinate subroutine
      TODO: The tstep conditional is appropriate if we assume there is no intervention during warmup
      It won't work if we introduce interventions into a scenario with a pre-existing intervention.
    */
    if ( Simulation::timeStep >  0) {
      nextDose=_lastVaccineDose+1;
      if (nextDose <= (int)Vaccine::_numberOfEpiDoses){
        if (W_UNIFORM() <  Vaccine::vaccineCoverage[nextDose - 1] && 
            Vaccine::targetagetstep[nextDose - 1] ==  agetstep) {
          vaccinate();
          Simulation::gMainSummary->reportEPIVaccination(ageGroup());
        }
      }
    }
  }
  if (IPTIntervention::IPT) {
    if ( Simulation::timeStep >  0) {
      // assumes 5-day intervals and Niakhar seasonality
      static int IPT_MIN_INTERVAL[9] = { 42, 48, 54, 60, 66, 36, 30, 24, 18 };
      static int IPT_MAX_INTERVAL[9] = { 60, 66, 72, 78, 82, 54, 48, 42, 42 };
      
      if (iptiEffect >= 14 && iptiEffect <= 22) {
        int yearInterval = (Global::modIntervalsPerYear(_simulationTime)-1);
        if (yearInterval >= IPT_MIN_INTERVAL[iptiEffect-14] &&
            yearInterval <  IPT_MAX_INTERVAL[iptiEffect-14])
          setLastSPDose();
      } else
        setLastSPDose();
    }
  }
}

void Human::setLastSPDose() {
  int agetstep=_simulationTime-_dateOfBirth;
  for (int i=1;i<=IPTIntervention::numberOfIPTiDoses; i++) {
    if (IPTIntervention::iptiTargetagetstep[i - 1] == agetstep) {
      double rnum=W_UNIFORM();
      if ( (rnum) <  IPTIntervention::iptiCoverage[i - 1]) {
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
  if (getInterventions().getIptiDescription().present() &&
      _simulationTime-_tLastTreatment >= 1 &&
      _simulationTime-_tLastTreatment <=  4) {
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
  Simulation::gMainSummary->addToPyrogenicThreshold(age, _morbidityModel->getPyrogenThres());
  Simulation::gMainSummary->addToSumX(age, log(_morbidityModel->getPyrogenThres()+1.0));
}

void initDrugAction () {
  //Part of this will be auto generated

}


int Human::ageGroup() const{
  return Simulation::gMainSummary->ageGroup(getAgeInYears());
}

double Human::getAgeInYears() const{
  return 1.0*((_simulationTime-_dateOfBirth)*Global::interval)/daysInYear;
}

double Human::getAgeInYears(int time) const{
  return 1.0*((time-_dateOfBirth)*Global::interval)/daysInYear;
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
  if (( agetstep*Global::interval >  20) && ( _simulationTime*Global::interval >  20)) {
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
  if (Global::modelVersion & CASE_MANAGEMENT_V2) {
    returnValue=true;
  }
  else {
    nextRegimen=_caseManagement->getNextRegimen(_simulationTime, entrypoint, _tLastTreatment, _latestRegimen);
    //doCM(1);	//FIXME: is this running as well as old case-management code?
    if (_caseManagement->getProbabilityGetsTreatment(nextRegimen-1)*_treatmentSeekingFactor > (W_UNIFORM())){
      if (Global::modelVersion & INCLUDES_PK_PD){
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
  if (Global::modelVersion & CASE_MANAGEMENT_V2) {
    returnValue=true;      
  }
  /*
    TODO: is there a reason we cannot first decide if the patient is treated, then increase latestRegimen
    and latestTreatment, instead of resetting it if not treated? (Can't think of one, but
    do we want to change this section of code rather than just introducing the new alternative (TS))
  */ 
  nextRegimen=_caseManagement->getNextRegimen(_simulationTime, Diagnosis::SEVERE_MALARIA, _tLastTreatment, _latestRegimen);
  //doCM(3);	//FIXME: is this running as well as old case-management code?
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
  
  if (q[8] < 1.0) {
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
  }
  
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
  else // assume true, so we don't get another else case (DH): if( q[8] >=  prandom)
  {
    _latestEvent.update(_simulationTime, agegroup, Diagnosis::SEVERE_MALARIA, Outcome::PARASITES_ARE_CLEARED_PATIENT_RECOVERS_INPATIENTS);
    _tLastTreatment = _simulationTime;
    _latestRegimen = nextRegimen;
    returnValue=true;      
    Simulation::gMainSummary->reportTreatment(agegroup,_latestRegimen);
  }
  return returnValue;
}

bool Human::defineEvent(){
  bool valdefineEvent =false;
  double ageYears=getAgeInYears();
  double prEpisode=_morbidityModel->getPEpisode(_timeStepMaxDensity,_totalDensity);
  
  //Fixed severe threshold
  double severeMalThreshold=sevMal_21+1;
  double prSevereEpisode=1-1/(1+_timeStepMaxDensity/severeMalThreshold);
  
  //Decide whether a clinical episode occurs and if so, which type
  double pCoinfection=comorbintercept_24/(1+ageYears/critAgeComorb_30);
  pCoinfection*=_comorbidityFactor;
  
  if ((W_UNIFORM())< prEpisode) {
    if (W_UNIFORM() < prSevereEpisode ||	// severe
        W_UNIFORM() < pCoinfection)		// coinfection
      valdefineEvent=severeMalaria();
    else
      valdefineEvent=uncomplicatedEvent(true);	// uncomplicated
    
    /*
      Indirect mortality	
      IndirectRisk is the probability of dying from indirect effects of malaria
      conditional on not having an acute attack of malaria
    */
    double indirectRisk=indirRiskCoFactor_18/(1+ageYears/critAgeComorb_30);
    indirectRisk=indirectRisk*_comorbidityFactor;
    if ( (W_UNIFORM()) < indirectRisk) {
      if (_doomed ==  0) {
        _doomed=-1;
      }
    }
    // Penalizing immunity for clinical episodes	
    if (Global::modelVersion & PENALISATION_EPISODES) {
      _cumulativeY=(double)_cumulativeYlag+(-(immPenalty_22*(_cumulativeY-_cumulativeYlag)));
      if (_cumulativeY <  0) {
        _cumulativeY=0.0;
      }
    }
  }
  else if(Global::modelVersion & NON_MALARIA_FEVERS) {
    double prNonMalariaFever=pCoinfection*RelativeRiskNonMalariaFever;
    if ((W_UNIFORM()) < prNonMalariaFever){
      valdefineEvent=uncomplicatedEvent(false);
    }
  }
  return valdefineEvent;
}


double Human::calculateSelectionCoefficient (Infection& inf) {
  return 1.0;
}

list<Drug*>* Human::getDrugs() {
  return &_drugs;
}

double Human::entoAvailability () {
  return _entoAvailability
      * entoInterventionITN.availability()
      * entoInterventionIRS.availability();
}
double Human::probMosqSurvivalBiting () {
  return _probMosqSurvivalBiting
      * entoInterventionITN.probMosqBiting();
}
double Human::probMosqSurvivalResting () {
  return _probMosqSurvivalResting
      * entoInterventionIRS.probMosqSurvivalResting();
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
  out << human._TBVEfficacy << endl; 
  out << human._totalDensity << endl; 
  out << human._ylag[0] << endl; 
  out << human._ylag[1] << endl; 
  out << human._ylag[2] << endl; 
  out << human._ylag[3] << endl; 
  out << human._lastSPDose << endl; 
  out << human._lastIptiOrPlacebo << endl; 
  out << human._comorbidityFactor << endl;  
  out << human._treatmentSeekingFactor << endl; 
  out << human._BaselineAvailabilityToMosquitoes << endl; 
  out << human._latestEvent;   
  out << *human._withinHostModel;
  human._morbidityModel->write(out);
  out << human._entoAvailability << endl;
  out << human._probMosqSurvivalBiting << endl;
  out << human._probMosqSurvivalResting << endl;
  out << human.entoInterventionITN;
  out << human.entoInterventionIRS;

  for(std::list<Infection*>::const_iterator iter=human.infections.begin(); iter != human.infections.end(); iter++)
    (*iter)->write (out);
     
  if (Global::modelVersion & INCLUDES_PK_PD) {
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
  in >> human._TBVEfficacy; 
  in >> human._totalDensity; 
  in >> human._ylag[0]; 
  in >> human._ylag[1]; 
  in >> human._ylag[2]; 
  in >> human._ylag[3]; 
  in >> human._lastSPDose; 
  in >> human._lastIptiOrPlacebo; 
  in >> human._comorbidityFactor;  
  in >> human._treatmentSeekingFactor; 
  in >> human._BaselineAvailabilityToMosquitoes; 
  in >> human._latestEvent;
  in >> *(human._withinHostModel);
  human._morbidityModel->read(in);
  in >> human._entoAvailability;
  in >> human._probMosqSurvivalBiting;
  in >> human._probMosqSurvivalResting;
  in >> human.entoInterventionITN;
  in >> human.entoInterventionIRS;

  for(int i=0;i<human._MOI;++i) {
    Infection* infection = new DescriptiveInfection();
    infection->read(in);
    human.infections.push_back(infection);
  }

  if (Global::modelVersion & INCLUDES_PK_PD) {
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
