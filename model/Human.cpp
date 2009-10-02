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
#include "Human.h"

#include <string>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <stdexcept>

#include "simulation.h"
#include "inputData.h"
#include "util/gsl.h"
#include "summary.h"
#include "intervention.h"
#include "Transmission/TransmissionModel.h"
#include "InfectionIncidenceModel.h"
#include "Clinical/ClinicalModel.h"
#include "WithinHost/DescriptiveIPT.h"	// only for summarizing


/*
  Constants common to all humans
  Decay in anti-toxic immunity
*/
//TODOConversion
#ifdef _WIN32
#define isnan(x) ((x) != (x))
#endif


// -----  Static functions  -----

void Human::initHumanParameters () {	// static
  // Init models used by humans:
  PerHostTransmission::initParameters(getInterventions());
  InfectionIncidenceModel::init();
  WithinHostModel::init();
  ClinicalModel::init();
  Vaccine::initParameters();
}

void Human::clear() {	// static clear
  WithinHostModel::clear();
  Vaccine::clearParameters();
}


// -----  Non-static functions: creation/destruction, checkpointing  -----

// Create new human
Human::Human(TransmissionModel& tm, int ID, int dateOfBirth, int simulationTime) :
    perHostTransmission(),
    infIncidence(InfectionIncidenceModel::createModel()),
    withinHostModel(WithinHostModel::createWithinHostModel()),
    _dateOfBirth(dateOfBirth), _ID(ID),
    _lastVaccineDose(0),
    _BSVEfficacy(0.0), _PEVEfficacy(0.0), _TBVEfficacy(0.0)
{
  if (_dateOfBirth != simulationTime && (Simulation::timeStep >= 0 || _dateOfBirth > simulationTime))
    throw out_of_range ("Invalid date of birth!");
  
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
    if (gsl::rngUniform() < 0.5) {
      availabilityFactor=1.8;
    }
  }
  if (Global::modelVersion & COMORB_HET) {
    _comorbidityFactor=0.2;
    if (gsl::rngUniform() < 0.5) {
      _comorbidityFactor=1.8;
    }	
  }
  if (Global::modelVersion & TREAT_HET) {
    _treatmentSeekingFactor=0.2;
    if (gsl::rngUniform() < 0.5) {            
      _treatmentSeekingFactor=1.8;
    }	
  }
  if (Global::modelVersion & TRANS_TREAT_HET) {
    _treatmentSeekingFactor=0.2;
    availabilityFactor=1.8;
    if (gsl::rngUniform()<0.5) {
      _treatmentSeekingFactor=1.8;
      availabilityFactor=0.2;
    }
  } else if (Global::modelVersion & COMORB_TRANS_HET) {
    if (gsl::rngUniform()<0.5) {
      _treatmentSeekingFactor=0.2;
    } else {
      _treatmentSeekingFactor=1.8;
    }
    availabilityFactor=1.8;
    _comorbidityFactor=1.8;
    if (gsl::rngUniform()<0.5) {
      availabilityFactor=0.2;
      _comorbidityFactor=0.2;
    }
  } else if (Global::modelVersion & TRIPLE_HET) {
    availabilityFactor=1.8;
    _comorbidityFactor=1.8;
    _treatmentSeekingFactor=0.2;
    if (gsl::rngUniform()<0.5) {
      availabilityFactor=0.2;
      _comorbidityFactor=0.2;
      _treatmentSeekingFactor=1.8;
    }
  }
  perHostTransmission.initialise (tm, availabilityFactor * infIncidence->getAvailabilityFactor(1.0));
  clinicalModel=ClinicalModel::createClinicalModel (_comorbidityFactor, _treatmentSeekingFactor);
}

// Load human from checkpoint
Human::Human(istream& in, TransmissionModel& tm) :
    perHostTransmission(in, tm),
    infIncidence(InfectionIncidenceModel::createModel(in)),
    withinHostModel(WithinHostModel::createWithinHostModel(in)),
    clinicalModel(ClinicalModel::createClinicalModel(in))
{
  in >> _dateOfBirth; 
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
  delete infIncidence;
  delete withinHostModel;
  delete clinicalModel;
}

void Human::write (ostream& out) const{
  perHostTransmission.write (out);
  infIncidence->write (out);
  withinHostModel->write (out);
  clinicalModel->write (out);
  out << _dateOfBirth << endl; 
  out << _ID << endl ; 
  out << _lastVaccineDose << endl;
  out << _BSVEfficacy << endl; 
  out << _PEVEfficacy << endl; 
  out << _TBVEfficacy << endl; 
  out << _ylag[0] << endl; 
  out << _ylag[1] << endl; 
  out << _ylag[2] << endl; 
  out << _ylag[3] << endl; 
}


// -----  Non-static functions: per-timestep update  -----

bool Human::update(int simulationTime, TransmissionModel* transmissionModel) {
  int ageTimeSteps = simulationTime-_dateOfBirth;
  if (clinicalModel->isDead(ageTimeSteps))
    return true;
  
  updateInterventionStatus(); 
  withinHostModel->updateImmuneStatus();
  updateInfection(transmissionModel);
  clinicalModel->update (*withinHostModel, getAgeInYears(), Simulation::simulationTime-_dateOfBirth);
  withinHostModel->update();
  clinicalModel->updateInfantDeaths (ageTimeSteps);
  _probTransmissionToMosquito = calcProbTransmissionToMosquito ();
  return false;
}

void Human::updateInfection(TransmissionModel* transmissionModel){
  int numInf = infIncidence->numNewInfections(transmissionModel->getEIR(Simulation::simulationTime, perHostTransmission, getAgeInYears()),
					      _PEVEfficacy, perHostTransmission);
  for (int i=1;i<=numInf; i++) {
    withinHostModel->newInfection();
  }
  
  withinHostModel->clearOldInfections();
  
  // _ylag is designed for a 5-day timestep model
  //TODO equiv for 1-day tstep
  if ((Simulation::simulationTime*Global::interval) % 5 == 0) {
    for (int i=3;i>0; i--) {
      _ylag[i]=_ylag[i-1];
    }
  }
  _ylag[0]=withinHostModel->getTotalDensity();
  
  withinHostModel->calculateDensities(getAgeInYears(), _BSVEfficacy);
}

void Human::updateInterventionStatus() {
  int ageTimeSteps = Simulation::simulationTime-_dateOfBirth;
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
	if (gsl::rngUniform() <  Vaccine::vaccineCoverage[_lastVaccineDose] &&
	    Vaccine::targetAgeTStep[_lastVaccineDose] == ageTimeSteps) {
          vaccinate();
          Simulation::gMainSummary->reportEPIVaccination(ageGroup());
        }
      }
    }
  }
  withinHostModel->IPTSetLastSPDose(ageTimeSteps, ageGroup());
  perHostTransmission.continousItnDistribution (ageTimeSteps);
}


void Human::massVaccinate () {
  vaccinate();
  Simulation::gMainSummary->reportMassVaccination(ageGroup());
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

void Human::IPTiTreatment () {
  withinHostModel->IPTiTreatment (ageGroup());
}


void Human::clearInfections () {
  //NOTE: if Population::massTreatment is incompatible with IPT, we can just pass false:
  withinHostModel->clearInfections(clinicalModel->latestDiagnosisIsSevereMalaria());
}

int Human::ageGroup() const{
  return Simulation::gMainSummary->ageGroup(getAgeInYears());
}

double Human::getAgeInYears() const{
  return double((Simulation::simulationTime-_dateOfBirth)*Global::interval) / daysInYear;
}


void Human::summarize(){
  if (DescriptiveIPTWithinHost::iptActive && clinicalModel->recentTreatment())
    return;	//NOTE: do we need this?
  
  double age = getAgeInYears();
  Simulation::gMainSummary->addToHost(age,1);
  withinHostModel->summarize(age);
  infIncidence->summarize (*Simulation::gMainSummary, age);
  clinicalModel->summarize (*Simulation::gMainSummary, age);
}


double Human::calcProbTransmissionToMosquito() const {
  /* This model is only really valid for 5-day timesteps; if it's needed for
  a one-day timestep, it should be redesigned; design of this model is
  described in AJTMH pp.32-33*/
  //FIXME: implement something valid for 1-day timestep
  // Also TODO: do we have to use a run-time switch here? Make interval compile-time constant?
  
  int ageTimeSteps=Simulation::simulationTime-_dateOfBirth;
  if (ageTimeSteps*Global::interval <= 20 || Simulation::simulationTime*Global::interval <= 20)
    return 0.0;
  
  //Infectiousness parameters: see AJTMH p.33, tau=1/sigmag**2 
  static const double beta1=1.0;
  static const double beta2=0.46;
  static const double beta3=0.17;
  static const double tau= 0.066;
  static const double mu= -8.1;
  
  // Take weighted sum of total asexual blood stage density 10, 15 and 20 days before:
  double x=beta1*_ylag[1]+beta2*_ylag[2]+beta3*_ylag[3];
  if (x < 0.001)
    return 0.0;
  
  double zval=(log(x)+mu)/sqrt(1.0/tau);
  double pone = gsl::cdfUGaussianP (zval);
  double transmit=(pone*pone);
  //transmit has to be between 0 and 1
  transmit=std::max(transmit, 0.0);
  transmit=std::min(transmit, 1.0);
  
  //	Include here the effect of transmission-blocking vaccination
  return transmit*(1.0-_TBVEfficacy);
}
