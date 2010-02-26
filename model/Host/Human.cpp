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
#include "Host/Human.h"

#include "Host/intervention.h"
#include "Host/InfectionIncidenceModel.h"
#include "Clinical/ClinicalModel.h"
#include "WithinHost/DescriptiveIPT.h"	// only for summarizing

#include "inputData.h"
#include "util/gsl.h"
#include "Transmission/TransmissionModel.h"
#include "Surveys.h"
#include "util/ModelOptions.hpp"

#include <string>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace OM { namespace Host {
    int Human::_ylagLen;


// -----  Static functions  -----

void Human::initHumanParameters () {	// static
  // Init models used by humans:
  Transmission::PerHostTransmission::initParameters(InputData.getInterventions());
  InfectionIncidenceModel::init();
  WithinHost::WithinHostModel::init();
  Clinical::ClinicalModel::init();
  Vaccine::initParameters();
  _ylagLen = Global::intervalsPer5Days * 4;
}

void Human::clear() {	// static clear
  WithinHost::WithinHostModel::clear();
  Vaccine::clearParameters();
}


// -----  Non-static functions: creation/destruction, checkpointing  -----

// Create new human
Human::Human(Transmission::TransmissionModel& tm, int dateOfBirth, int simulationTime) :
    perHostTransmission(),
    infIncidence(InfectionIncidenceModel::createModel()),
    withinHostModel(WithinHost::WithinHostModel::createWithinHostModel()),
    _dateOfBirth(dateOfBirth),
    _lastVaccineDose(0),
    _BSVEfficacy(0.0), _PEVEfficacy(0.0), _TBVEfficacy(0.0),
    _probTransmissionToMosquito(0.0)
{
  if (_dateOfBirth != simulationTime && (Global::simulationTime > 0 || _dateOfBirth > simulationTime))
    throw out_of_range ("Invalid date of birth!");
  
  _ylag.assign (_ylagLen, 0.0);
  
  
  /* Human heterogeneity; affects:
   * _comorbidityFactor (stored in PathogenesisModel)
   * _treatmentSeekingFactor (stored in CaseManagementModel)
   * availabilityFactor (stored in Transmission::PerHostTransmission)
   */
  double _comorbidityFactor = 1.0;
  double _treatmentSeekingFactor = 1.0;
  double availabilityFactor = 1.0;
  
  if (util::ModelOptions::option (util::TRANS_HET)) {
    availabilityFactor=0.2;
    if (rng::uniform01() < 0.5) {
      availabilityFactor=1.8;
    }
  }
  if (util::ModelOptions::option (util::COMORB_HET)) {
    _comorbidityFactor=0.2;
    if (rng::uniform01() < 0.5) {
      _comorbidityFactor=1.8;
    }	
  }
  if (util::ModelOptions::option (util::TREAT_HET)) {
    _treatmentSeekingFactor=0.2;
    if (rng::uniform01() < 0.5) {            
      _treatmentSeekingFactor=1.8;
    }	
  }
  if (util::ModelOptions::option (util::TRANS_TREAT_HET)) {
    _treatmentSeekingFactor=0.2;
    availabilityFactor=1.8;
    if (rng::uniform01()<0.5) {
      _treatmentSeekingFactor=1.8;
      availabilityFactor=0.2;
    }
  } else if (util::ModelOptions::option (util::COMORB_TRANS_HET)) {
    if (rng::uniform01()<0.5) {
      _treatmentSeekingFactor=0.2;
    } else {
      _treatmentSeekingFactor=1.8;
    }
    availabilityFactor=1.8;
    _comorbidityFactor=1.8;
    if (rng::uniform01()<0.5) {
      availabilityFactor=0.2;
      _comorbidityFactor=0.2;
    }
  } else if (util::ModelOptions::option (util::TRIPLE_HET)) {
    availabilityFactor=1.8;
    _comorbidityFactor=1.8;
    _treatmentSeekingFactor=0.2;
    if (rng::uniform01()<0.5) {
      availabilityFactor=0.2;
      _comorbidityFactor=0.2;
      _treatmentSeekingFactor=1.8;
    }
  }
  perHostTransmission.initialise (tm, availabilityFactor * infIncidence->getAvailabilityFactor(1.0));
  clinicalModel = Clinical::ClinicalModel::createClinicalModel (_comorbidityFactor, _treatmentSeekingFactor);
}

void Human::destroy() {
  delete infIncidence;
  delete withinHostModel;
  delete clinicalModel;
}


// -----  Non-static functions: per-timestep update  -----

bool Human::update(int simulationTime, Transmission::TransmissionModel* transmissionModel) {
  int ageTimeSteps = simulationTime-_dateOfBirth;
  if (clinicalModel->isDead(ageTimeSteps))
    return true;
  
  updateInterventionStatus();
  updateInfection(transmissionModel);
  clinicalModel->update (*withinHostModel, getAgeInYears(), Global::simulationTime-_dateOfBirth);
  clinicalModel->updateInfantDeaths (ageTimeSteps);
  _probTransmissionToMosquito = calcProbTransmissionToMosquito ();
  return false;
}

void Human::updateInfection(Transmission::TransmissionModel* transmissionModel){
  int numInf = infIncidence->numNewInfections(transmissionModel->getEIR(Global::simulationTime, perHostTransmission, getAgeInYears()),
					      _PEVEfficacy, perHostTransmission);
  for (int i=1;i<=numInf; i++) {
    withinHostModel->newInfection();
  }
  
  // Cache total density for infectiousness calculations
  _ylag[Global::simulationTime%_ylagLen]=withinHostModel->getTotalDensity();
  
  withinHostModel->calculateDensities(getAgeInYears(), _BSVEfficacy);
}

void Human::updateInterventionStatus() {
  int ageTimeSteps = Global::simulationTime-_dateOfBirth;
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
    
    //_ctsIntervs.deploy(ageTimeSteps);
    if (Global::timeStep >= 0) {
      if (_lastVaccineDose < (int)Vaccine::_numberOfEpiDoses){
	  if (Vaccine::targetAgeTStep[_lastVaccineDose] == ageTimeSteps &&
	      rng::uniform01() <  Vaccine::vaccineCoverage[_lastVaccineDose] ) {
          vaccinate();
          Surveys.current->reportEPIVaccinations (ageGroup(), 1);
        }
      }
    }
  }
  withinHostModel->IPTSetLastSPDose(ageTimeSteps, ageGroup());
  perHostTransmission.continousItnDistribution (ageTimeSteps);
}


void Human::massVaccinate () {
  vaccinate();
  Surveys.current->reportMassVaccinations (ageGroup(), 1);
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


void Human::massDrugAdministration () {
    clinicalModel->massDrugAdministration (*withinHostModel, getAgeInYears());
}

SurveyAgeGroup Human::ageGroup() const{
  return SurveyAgeGroup(getAgeInYears());
}

double Human::getAgeInYears() const{
  return double((Global::simulationTime-_dateOfBirth)*Global::interval) / Global::DAYS_IN_YEAR;
}


void Human::summarize(Survey& survey) {
  if (WithinHost::DescriptiveIPTWithinHost::iptActive && clinicalModel->recentTreatment())
    return;	//NOTE: this modifies the denominator to treat the 4*5 day intervals after an episode as 'not at risk' to match the IPTi trials
  
  SurveyAgeGroup ageGrp = ageGroup();
  survey.reportHosts (ageGrp, 1);
  withinHostModel->summarize (survey, ageGrp);
  infIncidence->summarize (survey, ageGrp);
  clinicalModel->summarize (survey, ageGrp);
}


double Human::calcProbTransmissionToMosquito() const {
  /* This model was designed for 5-day timesteps. We use the same model
  (sampling 10, 15 and 20 days ago) for 1-day timesteps to avoid having to
  (design and analyse a new model. Description: AJTMH pp.32-33 */
  int ageTimeSteps=Global::simulationTime-_dateOfBirth;
  if (ageTimeSteps*Global::interval <= 20 || Global::simulationTime*Global::interval <= 20)
    return 0.0;
  
  //Infectiousness parameters: see AJTMH p.33, tau=1/sigmag**2 
  static const double beta1=1.0;
  static const double beta2=0.46;
  static const double beta3=0.17;
  static const double tau= 0.066;
  static const double mu= -8.1;
  
  // Take weighted sum of total asexual blood stage density 10, 15 and 20 days before.
  // These values are one timestep more recent than that, however the calculated
  // value is not used until the next timestep when then ages would be correct.
  double x = beta1 * _ylag[(Global::simulationTime-2*Global::intervalsPer5Days+1) % _ylagLen]
	   + beta2 * _ylag[(Global::simulationTime-3*Global::intervalsPer5Days+1) % _ylagLen]
	   + beta3 * _ylag[(Global::simulationTime-4*Global::intervalsPer5Days+1) % _ylagLen];
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

} }