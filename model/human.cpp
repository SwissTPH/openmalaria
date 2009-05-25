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
#include "CaseManagementModel.h"
#include "summary.h"
#include "proteome.h"
#include "intervention.h"
#include "TransmissionModel.h"
#include <stdexcept>


/*
  Constants common to all humans
  Decay in anti-toxic immunity
*/
//TODOConversion
#ifdef _WIN32
#define isnan(x) ((x) != (x))
#endif


double Human::detectionlimit;

void Human::initHumanParameters () {	// static
  // Init models used by humans:
  WithinHostModel::init();
  PresentationModel::init();
  PerHostTransmission::initParameters();
  Vaccine::initParameters();
  
  /*
    Init parameters that are common to all humans
  */
  
  // TODO: Change values of these probabilities - should get from old C code for entomodel.
  // We assume that baseEntoAvailability is not the same as the availability used
  // by NonVector.cpp.

  double densitybias;
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
}

void Human::clear() {	// static clear
  WithinHostModel::clear();
  Vaccine::clearParameters();
}

// Create new human
Human::Human(TransmissionModel& tm, int ID, int dateOfBirth, int simulationTime) 
  : _perHostTransmission(tm), _simulationTime(simulationTime)
{
  //std::cout<<"newH:ID dateOfBirth "<<ID<<" "<<dateOfBirth<<std::endl;
  _BSVEfficacy=0.0;
  _dateOfBirth=dateOfBirth;
  if (_dateOfBirth > simulationTime) {
    // This test may be totally unnecessary; it was done in oldWithinHostModel
    throw out_of_range ("date of birth in future!");
  }
  _doomed=0;
  _timeStepMaxDensity=0.0;
  _ID=ID;
  //NOTE: initialized here to preserve order of random calls
  _withinHostModel = WithinHostModel::createWithinHostModel();
  _lastVaccineDose=0;
  _PEVEfficacy=0.0;
  _TBVEfficacy=0.0;
  _totalDensity=0.0;
  for (int i=1;i<=4; i++) {
    _ylag[i-1]=0.0;
  }
  
  // Stored in presentation model, not Human, now. Initialization looks somewhat entwined though..
  double _comorbidityFactor;
  // Ditto, but in case management.
  double _treatmentSeekingFactor;
  // ditto for _BaselineAvailabilityToMosquitoes in PerHostTransmission
  
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
    _perHostTransmission._BaselineAvailabilityToMosquitoes=1.8;
    if (W_UNIFORM()<0.5) {
      _treatmentSeekingFactor=1.8;
      _perHostTransmission._BaselineAvailabilityToMosquitoes=0.2;
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
    _perHostTransmission._BaselineAvailabilityToMosquitoes=1.8;
    _comorbidityFactor=1.8;
    if (W_UNIFORM()<0.5) {
      _perHostTransmission._BaselineAvailabilityToMosquitoes=0.2;
      _comorbidityFactor=0.2;
    }
  }  
  if (Global::modelVersion & TRIPLE_HET) {
    _perHostTransmission._BaselineAvailabilityToMosquitoes=1.8;
    _comorbidityFactor=1.8;
    _treatmentSeekingFactor=0.2;
    if (W_UNIFORM()<0.5) {
      _perHostTransmission._BaselineAvailabilityToMosquitoes=0.2;
      _comorbidityFactor=0.2;
      _treatmentSeekingFactor=1.8;
    }
  }
  _presentationModel=PresentationModel::createPresentationModel(_comorbidityFactor);
  _caseManagement = CaseManagementModel::createCaseManagementModel(_treatmentSeekingFactor);
}

// Load human from checkpoint
Human::Human(istream& in, TransmissionModel& tm, int simulationTime) 
  : _perHostTransmission(in, tm), _simulationTime(simulationTime)
{
  // NOTE: makes some unnecessary random calls
  // WARNING: this will likely change some tests with checkpointing
  _withinHostModel = WithinHostModel::createWithinHostModel();
  _presentationModel=PresentationModel::createPresentationModel(1.0);
  _caseManagement = CaseManagementModel::createCaseManagementModel(1.0);
  
  _withinHostModel->read(in);
  _presentationModel->read(in);
  _caseManagement->read (in);
  in >> _dateOfBirth; 
  in >> _doomed; 
  in >> _ID; 
  in >> _lastVaccineDose; 
  in >> _BSVEfficacy; 
  in >> _timeStepMaxDensity; 
  in >> _PEVEfficacy; 
  in >> _TBVEfficacy; 
  in >> _totalDensity; 
  in >> _ylag[0]; 
  in >> _ylag[1]; 
  in >> _ylag[2]; 
  in >> _ylag[3]; 
}

void Human::destroy() {
  delete _withinHostModel;
  delete _presentationModel;
  delete _caseManagement; 
}


void Human::updateInfection(double expectedInfectionRate, double expectedNumberOfInfections){
  int numInf = _perHostTransmission.numNewInfections(expectedInfectionRate, expectedNumberOfInfections);
  for (int i=1;i<=numInf; i++) {
    _withinHostModel->newInfection();
  }
  
  _withinHostModel->clearOldInfections();

  if ((_simulationTime*Global::interval) % 5 ==  0) {
    for (int i=1;i<=3; i++) {
      _ylag[4-i]=_ylag[3-i];
    }
  }
  _ylag[0]=_totalDensity;

  _withinHostModel->calculateDensities(*this);
}

bool Human::update(int simulationTime, TransmissionModel* transmissionModel) {
  int ageTimeSteps = simulationTime-getDateOfBirth();
  if (ageTimeSteps > Global::maxAgeIntervals) {	// too old
    _doomed = 1;
  }
  if (_doomed > 0)
    return true;	// remove from population
  
  _simulationTime = simulationTime;
  double expectedInfectionRate = transmissionModel->getRelativeAvailability(getAgeInYears()) * transmissionModel->getEIR(_simulationTime, _perHostTransmission);
  
  updateInterventionStatus(); 
  _withinHostModel->updateImmuneStatus();
  updateInfection(expectedInfectionRate,
                  transmissionModel->getExpectedNumberOfInfections(*this, expectedInfectionRate));
  determineClinicalStatus();
  _withinHostModel->update(getAgeInYears());
  
  // update array for the infant death rates
  if (ageTimeSteps <= (int)Global::intervalsPerYear){
    ++Global::infantIntervalsAtRisk[ageTimeSteps-1];
    if ((_doomed == 4) || (_doomed == -6) || (_doomed == 6)){
      ++Global::infantDeaths[ageTimeSteps-1];
    }
  }
  return false;
}

void Human::determineClinicalStatus(){ //TODO: this function should not do case management
  //	Countdown to indirect mortality
  if ( _doomed < 0) {
    _doomed--;
  }
  
  //indirect death: if this human's about to die, don't worry about further episodes:
  if (_doomed ==  -7) {	//clinical episode 6 intervals before
    _caseManagement->getEvent().update(_simulationTime, ageGroup(), Diagnosis::INDIRECT_MALARIA_DEATH, Outcome::INDIRECT_DEATH);
    /*
    doomed=7 is the code for indirect death, and 6 for neonatal death.
    Individuals with positive doomed values are removed at the start of
    the next time step. They cannot be removed immediately because 
    their deaths need to be counted.
    */
    _doomed  = 7;
    //Indirect Neonatal mortality
    return;
  }
  // Neonatal mortality:
  if(_simulationTime-_dateOfBirth == 1) {
    if (PresentationModel::eventNeonatalMortality()) {
      _caseManagement->getEvent().update(_simulationTime, ageGroup(), Diagnosis::INDIRECT_MALARIA_DEATH, Outcome::INDIRECT_DEATH);
      _doomed  = 6;
      return;
    }
  }
  
  /* infectionEvent determines whether there is an acute episode, or concomitant fever and
  then whether the episode is severe, uncomplicated or there is an indirect
  death.
  doCaseManagement clears infections if there was an effective treatment, or calls medicate,
  and decides whether the patient lives, has sequelae, or dies.
  */
  _caseManagement->doCaseManagement (_presentationModel->infectionEvent (getAgeInYears(), _totalDensity, _timeStepMaxDensity),
                                     *_withinHostModel,
                                     getAgeInYears(),
                                     _doomed);
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
  _withinHostModel->IPTSetLastSPDose(Simulation::simulationTime-_dateOfBirth, ageGroup());
}

void Human::clearInfections () {
  _withinHostModel->clearInfections(_caseManagement->getEvent());
}

void Human::IPTiTreatment (double compliance) {
  _withinHostModel->IPTiTreatment (compliance, ageGroup());
}


void Human::summarize(){
  double age = getAgeInYears();
  if (getInterventions().getIptiDescription().present() && _caseManagement->recentTreatment())
    return ;
  
  Simulation::gMainSummary->addToHost(age,1);
  _withinHostModel->summarize(age);
  if ( _totalDensity >  detectionlimit) {
    Simulation::gMainSummary->addToPatentHost(age, 1);
    Simulation::gMainSummary->addToSumLogDensity(age, log(_totalDensity));
  }
  _perHostTransmission.summarize (*Simulation::gMainSummary, age);
  Simulation::gMainSummary->addToPyrogenicThreshold(age, _presentationModel->getPyrogenThres());
  Simulation::gMainSummary->addToSumX(age, log(_presentationModel->getPyrogenThres()+1.0));
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


ostream& operator<<(ostream& out, const Human& human){
  human._perHostTransmission.write (out);
  human._withinHostModel->write(out);
  human._presentationModel->write(out);
  human._caseManagement->write (out);
  out << human._dateOfBirth << endl; 
  out << human._doomed << endl; 
  out << human._ID << endl ; 
  out << human._lastVaccineDose << endl;
  out << human._BSVEfficacy << endl; 
  out << human._timeStepMaxDensity << endl; 
  out << human._PEVEfficacy << endl; 
  out << human._TBVEfficacy << endl; 
  out << human._totalDensity << endl; 
  out << human._ylag[0] << endl; 
  out << human._ylag[1] << endl; 
  out << human._ylag[2] << endl; 
  out << human._ylag[3] << endl; 
  
  return out;
}
