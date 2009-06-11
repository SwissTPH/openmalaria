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
#include "human.h"

#include <string>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <stdexcept>

#include "simulation.h"
#include "inputData.h"
#include "GSLWrapper.h"
#include "summary.h"
#include "intervention.h"
#include "TransmissionModel.h"
#include "CaseManagementModel.h"
#include "InfectionIncidenceModel.h"


/*
  Constants common to all humans
  Decay in anti-toxic immunity
*/
//TODOConversion
#ifdef _WIN32
#define isnan(x) ((x) != (x))
#endif


void Human::initHumanParameters () {	// static
  // Init models used by humans:
  PerHostTransmission::initParameters();
  InfectionIncidenceModel::init();
  WithinHostModel::init();
  PathogenesisModel::init();
  Vaccine::initParameters();
}

void Human::clear() {	// static clear
  WithinHostModel::clear();
  Vaccine::clearParameters();
}

// Create new human
Human::Human(TransmissionModel& tm, int ID, int dateOfBirth, int simulationTime) :
    perHostTransmission(),
    infIncidence(InfectionIncidenceModel::createModel()),
    withinHostModel(WithinHostModel::createWithinHostModel())
{
  //std::cout<<"newH:ID dateOfBirth "<<ID<<" "<<dateOfBirth<<std::endl;
  _BSVEfficacy=0.0;
  _dateOfBirth=dateOfBirth;
  if (_dateOfBirth > simulationTime) {
    // This test may be totally unnecessary; it was done in oldWithinHostModel
    throw out_of_range ("date of birth in future!");
  }
  _doomed=0;
  _ID=ID;
  _lastVaccineDose=0;
  _PEVEfficacy=0.0;
  _TBVEfficacy=0.0;
  for (int i=0;i<4; i++) {
    _ylag[i]=0.0;
  }
  
  
  /* Human heterogeneity; affects:
   * _comorbidityFactor (stored in PathogenesisModel)
   * _treatmentSeekingFactor (stored in CaseManagementModel)
   * availabilityFactor (stored in PerHostTransmission)
   */
  double _comorbidityFactor = 1.0;
  double _treatmentSeekingFactor = 1.0;
  double availabilityFactor = 1.0;
  
  if (Global::modelVersion & TRANS_HET) {
    availabilityFactor=0.2;
    if (W_UNIFORM() < 0.5) {
      availabilityFactor=1.8;
    }
  }
  if (Global::modelVersion & COMORB_HET) {
    _comorbidityFactor=0.2;
    if (W_UNIFORM() < 0.5) {
      _comorbidityFactor=1.8;
    }	
  }
  if (Global::modelVersion & TREAT_HET) {
    _treatmentSeekingFactor=0.2;
    if (W_UNIFORM() < 0.5) {            
      _treatmentSeekingFactor=1.8;
    }	
  }
  if (Global::modelVersion & TRANS_TREAT_HET) {
    _treatmentSeekingFactor=0.2;
    availabilityFactor=1.8;
    if (W_UNIFORM()<0.5) {
      _treatmentSeekingFactor=1.8;
      availabilityFactor=0.2;
    }
  } else if (Global::modelVersion & COMORB_TRANS_HET) {
    if (W_UNIFORM()<0.5) {
      _treatmentSeekingFactor=0.2;
    } else {
      _treatmentSeekingFactor=1.8;
    }
    availabilityFactor=1.8;
    _comorbidityFactor=1.8;
    if (W_UNIFORM()<0.5) {
      availabilityFactor=0.2;
      _comorbidityFactor=0.2;
    }
  } else if (Global::modelVersion & TRIPLE_HET) {
    availabilityFactor=1.8;
    _comorbidityFactor=1.8;
    _treatmentSeekingFactor=0.2;
    if (W_UNIFORM()<0.5) {
      availabilityFactor=0.2;
      _comorbidityFactor=0.2;
      _treatmentSeekingFactor=1.8;
    }
  }
  perHostTransmission.initialise (tm, availabilityFactor * infIncidence->getAvailabilityFactor(1.0));
  pathogenesisModel=PathogenesisModel::createPathogenesisModel(_comorbidityFactor);
  caseManagement = CaseManagementModel::createCaseManagementModel(_treatmentSeekingFactor);
}

// Load human from checkpoint
Human::Human(istream& in, TransmissionModel& tm) :
    perHostTransmission(in, tm),
    infIncidence(InfectionIncidenceModel::createModel(in)),
    withinHostModel(WithinHostModel::createWithinHostModel(in)),
    pathogenesisModel(PathogenesisModel::createPathogenesisModel(in)),
    caseManagement(CaseManagementModel::createCaseManagementModel(in))
{
  in >> _dateOfBirth; 
  in >> _doomed; 
  in >> _ID; 
  in >> _lastVaccineDose; 
  in >> _BSVEfficacy; 
  in >> _PEVEfficacy; 
  in >> _TBVEfficacy; 
  in >> _ylag[0]; 
  in >> _ylag[1]; 
  in >> _ylag[2]; 
  in >> _ylag[3]; 
}

void Human::destroy() {
  delete withinHostModel;
  delete pathogenesisModel;
  delete caseManagement; 
}

ostream& operator<<(ostream& out, const Human& human){
  human.perHostTransmission.write (out);
  human.infIncidence->write (out);
  human.withinHostModel->write (out);
  human.pathogenesisModel->write (out);
  human.caseManagement->write (out);
  out << human._dateOfBirth << endl; 
  out << human._doomed << endl; 
  out << human._ID << endl ; 
  out << human._lastVaccineDose << endl;
  out << human._BSVEfficacy << endl; 
  out << human._PEVEfficacy << endl; 
  out << human._TBVEfficacy << endl; 
  out << human._ylag[0] << endl; 
  out << human._ylag[1] << endl; 
  out << human._ylag[2] << endl; 
  out << human._ylag[3] << endl; 
  
  return out;
}


void Human::updateInfection(TransmissionModel* transmissionModel){
  int numInf = infIncidence->numNewInfections(transmissionModel->getEIR(Simulation::simulationTime, perHostTransmission, getAgeInYears()),
					      _PEVEfficacy, perHostTransmission);
  for (int i=1;i<=numInf; i++) {
    withinHostModel->newInfection();
  }
  
  withinHostModel->clearOldInfections();
  
  // _ylag is designed for a 5-day timestep model
  if ((Simulation::simulationTime*Global::interval) % 5 == 0) {
    for (int i=3;i>0; i--) {
      _ylag[i]=_ylag[i-1];
    }
  }
  //NOTE: should this not also run only every 5 days?
  _ylag[0]=withinHostModel->getTotalDensity();
  
  withinHostModel->calculateDensities(*this);
}

bool Human::update(int simulationTime, TransmissionModel* transmissionModel) {
  int ageTimeSteps = simulationTime-_dateOfBirth;
  if (ageTimeSteps > Global::maxAgeIntervals) {	// too old
    _doomed = 1;
  }
  if (_doomed > 0)
    return true;	// remove from population
  
  updateInterventionStatus(); 
  withinHostModel->updateImmuneStatus();
  updateInfection(transmissionModel);
  determineClinicalStatus();
  withinHostModel->update(getAgeInYears());
  
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
    caseManagement->getEvent().update(Simulation::simulationTime, ageGroup(), Diagnosis::INDIRECT_MALARIA_DEATH, Outcome::INDIRECT_DEATH);
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
  if(Simulation::simulationTime-_dateOfBirth == 1) {
    if (PathogenesisModel::eventNeonatalMortality()) {
      caseManagement->getEvent().update(Simulation::simulationTime, ageGroup(), Diagnosis::INDIRECT_MALARIA_DEATH, Outcome::INDIRECT_DEATH);
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
  caseManagement->doCaseManagement (pathogenesisModel->infectionEvent (getAgeInYears(), withinHostModel->getTotalDensity(), withinHostModel->getTimeStepMaxDensity()),
                                     *withinHostModel,
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
  if (Vaccine::anyVaccine) {
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
    if (Simulation::timeStep >= 0) {
      if (_lastVaccineDose < (int)Vaccine::_numberOfEpiDoses){
	if (W_UNIFORM() <  Vaccine::vaccineCoverage[_lastVaccineDose] &&
            Vaccine::targetagetstep[_lastVaccineDose] == Simulation::simulationTime-_dateOfBirth) {
          vaccinate();
          Simulation::gMainSummary->reportEPIVaccination(ageGroup());
        }
      }
    }
  }
  withinHostModel->IPTSetLastSPDose(Simulation::simulationTime-_dateOfBirth, ageGroup());
}

void Human::clearInfections () {
  withinHostModel->clearInfections(caseManagement->getEvent());
}

void Human::IPTiTreatment (double compliance) {
  withinHostModel->IPTiTreatment (compliance, ageGroup());
}


void Human::summarize(){
  double age = getAgeInYears();
  if (getInterventions().getIptiDescription().present() && caseManagement->recentTreatment())
    return ;
  
  Simulation::gMainSummary->addToHost(age,1);
  withinHostModel->summarize(age);
  infIncidence->summarize (*Simulation::gMainSummary, age);
  Simulation::gMainSummary->addToPyrogenicThreshold(age, pathogenesisModel->getPyrogenThres());
  Simulation::gMainSummary->addToSumX(age, log(pathogenesisModel->getPyrogenThres()+1.0));
}


int Human::ageGroup() const{
  return Simulation::gMainSummary->ageGroup(getAgeInYears());
}

double Human::getAgeInYears() const{
  return 1.0*((Simulation::simulationTime-_dateOfBirth)*Global::interval)/daysInYear;
}


double Human::infectiousness(){
  double transmit;
  //Infectiousness parameters: see AJTMH p.33, tau=1/sigmag**2 
  static const double beta2=0.46;
  static const double beta3=0.17;
  static const double tau= 0.066;
  static const double mu= -8.1;
  int agetstep=Simulation::simulationTime-_dateOfBirth;
  /*
    Original infectiousness model based on 5 day intervals updates
    lagged variables only every 5 days and cannot compute infectiousness
    for the first 20 days of the simulation
  */
  if ((agetstep*Global::interval >  20) && (Simulation::simulationTime*Global::interval >  20)) {
    double x=_ylag[1]+beta2*_ylag[2]+beta3*_ylag[3];
    if ( x <  0.001) {
      transmit=0.0;
    }
    else {
      double zval=(log(x)+mu)/sqrt(1.0/tau);
      double pone=W_UGAUSS_P(zval);
      transmit=(pone*pone);
      //transmit has to be between 0 and 1
      transmit=std::max(transmit, 0.0);
      transmit=std::min(transmit, 1.0);
    }
  } else {
    transmit=0.0;
  }
  //	Include here the effect of transmission-blocking vaccination
  return transmit*(1.0-_TBVEfficacy);
}
