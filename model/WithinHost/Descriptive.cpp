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

#include "util/gsl.h"
#include "human.h"
#include "WithinHost/Descriptive.h"
#include "simulation.h"
#include "intervention.h"
#include "summary.h"

using namespace std;

const int DescriptiveWithinHostModel::MAX_INFECTIONS = 21;


// -----  Initialization  -----

DescriptiveWithinHostModel::DescriptiveWithinHostModel() :
    WithinHostModel(), _MOI(0),
    _cumulativeY(0.0), _cumulativeh(0.0), _cumulativeYlag(0.0),
    patentInfections(0)
{
  _innateImmunity = gsl::rngGauss(0, sigma_i);
}

DescriptiveWithinHostModel::DescriptiveWithinHostModel(istream& in) :
    WithinHostModel(in)
{
  readDescriptiveWHM (in);
  
  for(int i=0;i<_MOI;++i)
    infections.push_back(new DescriptiveInfection(in));
}

DescriptiveWithinHostModel::DescriptiveWithinHostModel(istream& in, bool) :
    WithinHostModel(in)
{
  readDescriptiveWHM (in);
}

DescriptiveWithinHostModel::~DescriptiveWithinHostModel() {
  clearAllInfections();
}


// -----  Data checkpointing  -----

void DescriptiveWithinHostModel::write(ostream& out) const {
  writeDescriptiveWHM (out);
}

void DescriptiveWithinHostModel::readDescriptiveWHM (istream& in) {
  in >> _MOI;
  in >> patentInfections; 
  in >> _cumulativeh;
  in >> _cumulativeY;
  in >> _cumulativeYlag;
  in >> _innateImmunity;
  
  if (_MOI < 0 || _MOI > MAX_INFECTIONS)
    throw checkpoint_error ("_MOI");
}

void DescriptiveWithinHostModel::writeDescriptiveWHM(ostream& out) const {
  out << _cumulativeInfections << endl;
  out << _pTransToMosq << endl;
  out << totalDensity << endl;
  out << timeStepMaxDensity << endl;
  
  out << _MOI << endl;
  out << patentInfections << endl;
  out << _cumulativeh << endl;
  out << _cumulativeY << endl;
  out << _cumulativeYlag << endl;
  out << _innateImmunity << endl;
  
  for(std::list<DescriptiveInfection*>::const_iterator iter=infections.begin(); iter != infections.end(); iter++)
    (*iter)->write (out);
}

// -----  Update function, called each step  -----

void DescriptiveWithinHostModel::update () {}


// -----  Simple infection adders/removers  -----

void DescriptiveWithinHostModel::newInfection(){
  if (_MOI < MAX_INFECTIONS) {
    _cumulativeInfections++;
    infections.push_back(new DescriptiveInfection(Simulation::simulationTime));
    _MOI++;
  }
}

void DescriptiveWithinHostModel::clearOldInfections(){
  std::list<DescriptiveInfection*>::iterator iter=infections.begin();
  while(iter != infections.end()){
    if (Simulation::simulationTime >= (*iter)->getEndDate()) {
      delete *iter;
      iter=infections.erase(iter);
      _MOI--;
    }
    else{
      iter++;
    }
  }
}

void DescriptiveWithinHostModel::clearAllInfections(){
  std::list<DescriptiveInfection*>::iterator i;
  for(i=infections.begin(); i != infections.end(); i++){
    delete *i;
  }
  infections.clear();
  _MOI=0;
}

void DescriptiveWithinHostModel::medicate(string drugName, double qty, int time, double age) {}


// -----  immunity  -----

void DescriptiveWithinHostModel::updateImmuneStatus(){
  if (immEffectorRemain < 1){
    _cumulativeh*=immEffectorRemain;
    _cumulativeY*=immEffectorRemain;
  }
  if (asexImmRemain < 1){
    _cumulativeh*=asexImmRemain/
        (1+(_cumulativeh*(1-asexImmRemain)/Infection::cumulativeHstar));
    _cumulativeY*=asexImmRemain/
        (1+(_cumulativeY*(1-asexImmRemain)/Infection::cumulativeYstar));
  }
}

void DescriptiveWithinHostModel::immunityPenalisation() {
  _cumulativeY=(double)_cumulativeYlag-(immPenalty_22*(_cumulativeY-_cumulativeYlag));
  if (_cumulativeY <  0) {
    _cumulativeY=0.0;
  }
}


// -----  Density calculations  -----

// NOTE: refering back to human so much isn't good programming practice. Could
// some variables be stored locally?
void DescriptiveWithinHostModel::calculateDensities(Human& human) {
  _cumulativeYlag = _cumulativeY;
  
  double ageyears = human.getAgeInYears();
  _pTransToMosq = 0.0;
  patentInfections = 0;
  totalDensity = 0.0;
  timeStepMaxDensity = 0.0;
  if (_cumulativeInfections >  0) {
    // Values of _cumulativeh/Y at beginning of step
    // (values are adjusted for each infection)
    double cumulativeh=_cumulativeh;
    double cumulativeY=_cumulativeY;
    
    // IPTi SP dose clears infections at the time that blood-stage parasites appear     
    SPAction(human);
    
    std::list<DescriptiveInfection*>::iterator iter;
    for(iter=infections.begin(); iter!=infections.end(); iter++){
      //std::cout<<"uis: "<<infData->duration<<std::endl;
      double infStepMaxDens = timeStepMaxDensity;
      
      if (Global::modelVersion & MAX_DENS_RESET) {
        infStepMaxDens=0.0;
      }
      (*iter)->determineDensities(Simulation::simulationTime, cumulativeY, ageyears, cumulativeh , infStepMaxDens);
      (*iter)->multiplyDensity(exp(-_innateImmunity));

        /*
      Possibly a better model version ensuring that the effect of variation in innate immunity
          is reflected in case incidence would have the following here:
        */
      if (Global::modelVersion & INNATE_MAX_DENS) {
        infStepMaxDens *= exp(-_innateImmunity);
      }
        //Include here the effect of blood stage vaccination
      if (Vaccine::BSV.active) {
	double factor = 1.0-human.getBSVEfficacy();
	(*iter)->multiplyDensity(factor);
	infStepMaxDens *= factor;
      }
      
      // Include here the effect of attenuated infections by SP concentrations
      IPTattenuateAsexualDensity (**iter);
      
      if (Global::modelVersion & MAX_DENS_CORRECTION) {
        infStepMaxDens = std::max(infStepMaxDens, timeStepMaxDensity);
      }
      timeStepMaxDensity = infStepMaxDens;
      
      totalDensity += (*iter)->getDensity();
      //Compute the proportion of parasites remaining after innate blood stage effect
      if ((*iter)->getDensity() > detectionLimit) {
        patentInfections++;
      }
      if ((*iter)->getStartDate() == (Simulation::simulationTime-1)) {
        _cumulativeh++;
      }
      (*iter)->setDensity(std::min(maxDens, (*iter)->getDensity()));
      (*iter)->setCumulativeExposureJ((*iter)->getCumulativeExposureJ()+Global::interval*(*iter)->getDensity());
      _cumulativeY += Global::interval*(*iter)->getDensity();
    }
    
    IPTattenuateAsexualMinTotalDensity(human);
  }
  _pTransToMosq = human.infectiousness();
}

void DescriptiveWithinHostModel::SPAction(Human&){}
void DescriptiveWithinHostModel::IPTattenuateAsexualDensity (DescriptiveInfection&) {}
void DescriptiveWithinHostModel::IPTattenuateAsexualMinTotalDensity (Human&) {}

// -----  Summarize  -----

void DescriptiveWithinHostModel::summarize(double age) {
  if (_MOI > 0) {
    Simulation::gMainSummary->addToInfectedHost(age,1);
    Simulation::gMainSummary->addToTotalInfections(age, _MOI);
    Simulation::gMainSummary->addToTotalPatentInfections(age, patentInfections);
  }
  if (parasiteDensityDetectible()) {
    Simulation::gMainSummary->addToPatentHost(age, 1);
    Simulation::gMainSummary->addToSumLogDensity(age, log(totalDensity));
  }
}
