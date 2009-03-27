/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "OldIPTWithinHostModel.h"
#include "human.h"
#include "simulation.h"
#include "GSLWrapper.h"
#include "event.h"
#include "summary.h"
#include "intervention.h"


// -----  init  -----
OldIPTWithinHostModel::OldIPTWithinHostModel () :
    _SPattenuationt(0), _lastSPDose (missing_value), _lastIptiOrPlacebo (missing_value)
{
  if (Global::modelVersion & INCLUDES_PK_PD) {
    cerr << "OldIPTWithinHostModel not intended to work with DrugAction" << endl;
    // The IPT code has its own implementation of non-instantaneous drug action (SPAction, etc).
    throw 0;
  }
}


// -----  Simple infection adders/removers  -----

void OldIPTWithinHostModel::newInfection(){
  //std::cout<<"MOI "<<_MOI<<std::endl;
  if (_MOI <=  20) {
    _cumulativeInfections++;
    infections.push_back(DescriptiveInfection(_lastSPDose, Simulation::simulationTime));
    _MOI++;
  }
}


// -----    -----

void OldIPTWithinHostModel::IPTClearInfections (Event& latestEvent) {
  int fortnight = int((14.0/Global::interval)+0.5);	// round to nearest
  if ( latestEvent.getDiagnosis() ==  Diagnosis::SEVERE_MALARIA) {
    clearAllInfections();
  }
  else if(Simulation::simulationTime-_lastIptiOrPlacebo <=  fortnight) {
    clearAllInfections();
          // IPTi trials used quinine for fevers within 14 days of an ipti or placebo dose   
  }
  else if(Simulation::simulationTime-_lastSPDose <=  fortnight) {
    clearAllInfections();
          /*
    second line used if fever within 14 days of SP dose (ipti or treatment)
    TODO: if this code is to survive, then the iptiEffect values should be 
    symbolic constants
          */
  }

  else if( IPTIntervention::iptiEffect ==  2 ||  IPTIntervention::iptiEffect ==  12) {
    clearAllInfections();
    _lastSPDose=Simulation::simulationTime+1;
  }
  else if( IPTIntervention::iptiEffect ==  3 ||  IPTIntervention::iptiEffect ==  13) {
    clearAllInfections();
  }
  else if(IPTIntervention::iptiEffect >=  14 && IPTIntervention::iptiEffect < 30) {
    clearAllInfections();
  }
  else {
    _lastSPDose=Simulation::simulationTime+1;
    clearAllInfections();
          // SPAction will first act at the beginning of the next Global::interval
  }
}

void OldIPTWithinHostModel::IPTSetLastSPDose (int agetstep, int ageGroup) {
  if (Simulation::timeStep <= 0) return;
  
  // assumes 5-day intervals and Niakhar seasonality
  static int IPT_MIN_INTERVAL[9] = { 42, 48, 54, 60, 66, 36, 30, 24, 18 };
  static int IPT_MAX_INTERVAL[9] = { 60, 66, 72, 78, 82, 54, 48, 42, 42 };
  
  if (IPTIntervention::iptiEffect >= 14 && IPTIntervention::iptiEffect <= 22) {
    int yearInterval = (Global::modIntervalsPerYear(Simulation::simulationTime)-1);
    if (yearInterval <  IPT_MIN_INTERVAL[IPTIntervention::iptiEffect-14] &&
        yearInterval >= IPT_MAX_INTERVAL[IPTIntervention::iptiEffect-14])
      return;
  }
  
  for (int i=0;i<IPTIntervention::numberOfIPTiDoses; i++) {
    if (IPTIntervention::iptiTargetagetstep[i] == agetstep) {
      if (W_UNIFORM() <  IPTIntervention::iptiCoverage[i]) {
        _lastIptiOrPlacebo=Simulation::simulationTime;
        /*
        iptiEffect denotes treatment or placebo group
        and also the treatment given when sick (trial-dependent)
        */
        if (IPTIntervention::iptiEffect >=  10) {
          _lastSPDose=Simulation::simulationTime;
          Simulation::gMainSummary->reportIPTDose(ageGroup);
        }
      }
    }
  }
}

void OldIPTWithinHostModel::IPTiTreatment (double compliance, int ageGroup) {
  //NOTE to AR remove _cumulativeInfections > 0 &&
  if (_cumulativeInfections > 0 && W_UNIFORM() < compliance){
    _lastIptiOrPlacebo = Simulation::simulationTime;
    /*
    * iptiEffect denotes treatment or placebo group
    * and also the treatment given when sick (trial-dependent)
    */
    if (IPTIntervention::iptiEffect >= 10){
      _lastSPDose = Simulation::simulationTime;
      Simulation::gMainSummary->reportIPTDose(ageGroup);
    }
  }
}


// -----  density calculation  -----

void OldIPTWithinHostModel::SPAction(Human& human){
  /*TODO if we want to look at presumptive SP treatment with the PkPD model we
  need to add some code here that will be conditionally implemented depending on the
  model version.*/

  double rnum;
  std::list<DescriptiveInfection>::iterator i=infections.begin();
  while(i != infections.end()){
    if ( 1+Simulation::simulationTime-i->getStartDate()-Global::latentp > 0){
      rnum=W_UNIFORM();
      if ((rnum<=IPTIntervention::genotypeACR[i->getGenoTypeID()-1]) &&
           (Simulation::simulationTime - _lastSPDose <= IPTIntervention::genotypeProph[i->getGenoTypeID()-1])) {
        i->destroy();
        i=infections.erase(i);
        _MOI--;
           }
           else{
             i++;
           }
    }
    else{
      i++;
    }
  }
}

void OldIPTWithinHostModel::IPTattenuateAsexualDensity (std::list<DescriptiveInfection>::iterator i) {
  if (Global::modelVersion & ATTENUATION_ASEXUAL_DENSITY) {
    if (i->getSPattenuate() == 1) {
      i->multiplyDensity(1.0/IPTIntervention::genotypeAtten[i->getGenoTypeID() - 1]);
      timeStepMaxDensity=(double)timeStepMaxDensity/IPTIntervention::genotypeAtten[i->getGenoTypeID() - 1];
      _SPattenuationt=(int)std::max(_SPattenuationt*1.0, (i->getStartDate()+(i->getDuration()/Global::interval) * IPTIntervention::genotypeAtten[i->getGenoTypeID() - 1]));
    }
  }
}

void OldIPTWithinHostModel::IPTattenuateAsexualMinTotalDensity (Human& human) {
  if (Global::modelVersion & ATTENUATION_ASEXUAL_DENSITY) {
    if (_SPattenuationt > Simulation::simulationTime &&  human.getTotalDensity() <  10) {
      human.setTotalDensity(10);
      _cumulativeY += 10;
    }
  }
}


// -----  Data checkpointing  -----

void OldIPTWithinHostModel::read(istream& in) {
  readOWHM (in);
  in >> _SPattenuationt;
  in >> _lastSPDose; 
  in >> _lastIptiOrPlacebo; 
}
void OldIPTWithinHostModel::write(ostream& out) const {
  writeOWHM (out);
  out << _SPattenuationt << endl;
  out << _lastSPDose << endl; 
  out << _lastIptiOrPlacebo << endl; 
}
