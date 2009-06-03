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

#include "WithinHostModel/OldIPT.h"
#include "human.h"
#include "simulation.h"
#include "GSLWrapper.h"
#include "event.h"
#include "summary.h"
#include "intervention.h"
#include "inputData.h"
#include "WithinHostModel/OldIPTInfection.h"

bool OldIPTWithinHostModel::iptActive = false;

// -----  static data  -----

int OldIPTWithinHostModel::numberOfIPTiDoses = 0;
int *OldIPTWithinHostModel::iptiTargetagetstep;
double *OldIPTWithinHostModel::iptiCoverage;
int OldIPTWithinHostModel::iptiEffect;


// -----  init  -----

void OldIPTWithinHostModel::initParameters () {
  const scnXml::Interventions& xmlInterventions = getInterventions();
  iptActive = xmlInterventions.getIptiDescription().present();
  if (!iptActive) return;
  
  if (Global::interval != 5)
    throw domain_error ("IPT code only supports using an interval of 5");
  
  // --- IptiDescription begin ---
  const scnXml::IptDescription& xmlIPTI = xmlInterventions.getIptiDescription().get();
  
  iptiEffect = xmlIPTI.getIptiEffect();
  // --- IptiDescription end ---
  
  if (xmlInterventions.getContinuous().present()) {
    const scnXml::Continuous::IptiSequence& xmlIpti = xmlInterventions.getContinuous().get().getIpti();
    numberOfIPTiDoses = xmlIpti.size();
    
    iptiTargetagetstep = (int*)malloc(((numberOfIPTiDoses))*sizeof(int));
    iptiCoverage = (double*)malloc(((numberOfIPTiDoses))*sizeof(double));
    for (int i=0;i<numberOfIPTiDoses; i++) {
      iptiTargetagetstep[i] = floor(xmlIpti[i].getTargetAgeYrs() * daysInYear / (1.0*Global::interval));
      iptiCoverage[i] = xmlIpti[i].getCoverage();
    }
  } else
    numberOfIPTiDoses = 0;
  
  OldIPTInfection::initParameters(xmlInterventions);
}

void OldIPTWithinHostModel::clearParameters () {
  if (!iptActive) return;
  if (numberOfIPTiDoses) {
    free(iptiTargetagetstep);
    free(iptiCoverage);
  }
  OldIPTInfection::clearParameters();
}

OldIPTWithinHostModel::OldIPTWithinHostModel () :
    _SPattenuationt(TIMESTEP_NEVER), _lastSPDose (TIMESTEP_NEVER), _lastIptiOrPlacebo (TIMESTEP_NEVER)
{
  if (Global::modelVersion & INCLUDES_PK_PD) {
    throw xml_scenario_error ("OldIPTWithinHostModel not intended to work with DrugAction");
    // The IPT code has its own implementation of non-instantaneous drug action (SPAction, etc).
  }
}


// -----  Data checkpointing  -----

OldIPTWithinHostModel::OldIPTWithinHostModel (istream& in) :
    DescriptiveWithinHostModel(in, true)
{
  for(int i=0;i<_MOI;++i)
    infections.push_back(new OldIPTInfection(in));
  
  in >> _SPattenuationt;
  in >> _lastSPDose; 
  in >> _lastIptiOrPlacebo; 
}

void OldIPTWithinHostModel::write(ostream& out) const {
  writeDescriptiveWHM (out);
  
  out << _SPattenuationt << endl;
  out << _lastSPDose << endl; 
  out << _lastIptiOrPlacebo << endl; 
}


// -----  Simple infection adders/removers  -----

void OldIPTWithinHostModel::newInfection(){
  if (_MOI <= MAX_INFECTIONS) {
    _cumulativeInfections++;
    infections.push_back(new OldIPTInfection(_lastSPDose, Simulation::simulationTime));
    _MOI++;
  }
}

// -----    -----

void OldIPTWithinHostModel::clearInfections (Event& latestEvent) {
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

  else if( iptiEffect ==  2 ||  iptiEffect ==  12) {
    clearAllInfections();
    _lastSPDose=Simulation::simulationTime+1;
  }
  else if( iptiEffect ==  3 ||  iptiEffect ==  13) {
    clearAllInfections();
  }
  else if(iptiEffect >=  14 && iptiEffect < 30) {
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
  // These numbers, should have MAX = MIN + 18 (modulo 73).
  static int IPT_MIN_INTERVAL[9] = { 43, 49, 55, 61, 67, 37, 31, 25, 19 };
  static int IPT_MAX_INTERVAL[9] = { 61, 67, 73,  6, 12, 55, 49, 43, 31 };
  
  if (iptiEffect >= 14 && iptiEffect <= 22) {
    int yearInterval = (Global::modIntervalsPerYear(Simulation::simulationTime));
    //int yearInterval = Simulation::simulationTime % Global::intervalsPerYear;
    int min = IPT_MIN_INTERVAL[iptiEffect-14];
    int max = IPT_MAX_INTERVAL[iptiEffect-14];
    // We're using modular arithmatic here, to represent a time period 5*18 days long.
    if (min < max) {
      if (yearInterval < min || yearInterval >= max)
	return;
    } else {
      if (yearInterval < min && yearInterval >= max)
	return;
    }
  }
  
  for (int i=0;i<numberOfIPTiDoses; i++) {
    if (iptiTargetagetstep[i] == agetstep) {
      if (W_UNIFORM() <  iptiCoverage[i]) {
        _lastIptiOrPlacebo=Simulation::simulationTime;
        /*
        iptiEffect denotes treatment or placebo group
        and also the treatment given when sick (trial-dependent)
        */
        if (iptiEffect >=  10) {
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
    if (iptiEffect >= 10){
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

  std::list<DescriptiveInfection*>::iterator iter=infections.begin();
  while(iter != infections.end()) {
    if (1 + Simulation::simulationTime - (*iter)->getStartDate() > Global::latentp) {
      OldIPTInfection* infec = dynamic_cast<OldIPTInfection*> (*iter);
      if (infec == NULL) throw logic_error ("infections should be of type DescriptiveInfection");
      size_t genoTypeIndex = infec->getGenoTypeID()-1;
      if ((W_UNIFORM() <= OldIPTInfection::genotypeACR[genoTypeIndex]) &&
	  (Simulation::simulationTime - _lastSPDose <= OldIPTInfection::genotypeProph[genoTypeIndex])) {
        delete *iter;
        iter=infections.erase(iter);
        _MOI--;
      } else {
	iter++;
      }
    } else {
      iter++;
    }
  }
}

void OldIPTWithinHostModel::IPTattenuateAsexualDensity (DescriptiveInfection& dinfec) {
  if (!(Global::modelVersion & ATTENUATION_ASEXUAL_DENSITY)) return;
  
  OldIPTInfection& infec = dynamic_cast<OldIPTInfection&> (dinfec);
    if (infec.getSPattenuate() == 1) {
      infec.multiplyDensity(1.0/OldIPTInfection::genotypeAtten[infec.getGenoTypeID() - 1]);
      timeStepMaxDensity=(double)timeStepMaxDensity/OldIPTInfection::genotypeAtten[infec.getGenoTypeID() - 1];
      _SPattenuationt=(int)std::max(_SPattenuationt*1.0, (infec.getStartDate()+(infec.getDuration()/Global::interval) * OldIPTInfection::genotypeAtten[infec.getGenoTypeID() - 1]));
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
