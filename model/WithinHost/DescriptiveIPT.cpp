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

#include "WithinHost/DescriptiveIPT.h"
#include "WithinHost/DescriptiveIPTInfection.h"
#include "Simulation.h"
#include "util/gsl.h"
#include "summary.h"
#include "intervention.h"
#include "inputData.h"

bool DescriptiveIPTWithinHost::iptActive = false;

// -----  static data  -----

int DescriptiveIPTWithinHost::numberOfIPTiDoses = 0;
int *DescriptiveIPTWithinHost::iptiTargetagetstep;
double *DescriptiveIPTWithinHost::iptiCoverage;
int DescriptiveIPTWithinHost::iptiEffect;


// -----  init  -----

void DescriptiveIPTWithinHost::initParameters () {
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
      iptiTargetagetstep[i] = static_cast<int>(floor(xmlIpti[i].getTargetAgeYrs() * daysInYear / (1.0*Global::interval)));
      iptiCoverage[i] = xmlIpti[i].getCoverage();
    }
  } else
    numberOfIPTiDoses = 0;
  
  DescriptiveIPTInfection::initParameters(xmlInterventions);
}

void DescriptiveIPTWithinHost::clearParameters () {
  if (!iptActive) return;
  if (numberOfIPTiDoses) {
    free(iptiTargetagetstep);
    free(iptiCoverage);
  }
  DescriptiveIPTInfection::clearParameters();
}

DescriptiveIPTWithinHost::DescriptiveIPTWithinHost () :
    _SPattenuationt(TIMESTEP_NEVER), _lastSPDose (TIMESTEP_NEVER), _lastIptiOrPlacebo (TIMESTEP_NEVER),
    _cumulativeInfections(0)
{
  if (Global::modelVersion & INCLUDES_PK_PD) {
    throw xml_scenario_error ("DescriptiveIPTWithinHost not intended to work with DrugAction");
    // The IPT code has its own implementation of non-instantaneous drug action (SPAction, etc).
  }
}


// -----  Data checkpointing  -----

DescriptiveIPTWithinHost::DescriptiveIPTWithinHost (istream& in) :
    DescriptiveWithinHostModel(in, true)
{
  for(int i=0;i<_MOI;++i)
    infections.push_back(new DescriptiveIPTInfection(in));
  
  in >> _SPattenuationt;
  in >> _lastSPDose; 
  in >> _lastIptiOrPlacebo;
  in >> _cumulativeInfections;
}

void DescriptiveIPTWithinHost::write(ostream& out) const {
  DescriptiveWithinHostModel::write (out);
  
  out << _SPattenuationt << endl;
  out << _lastSPDose << endl; 
  out << _lastIptiOrPlacebo << endl;
  out << _cumulativeInfections << endl;
}


// -----  Simple infection adders/removers  -----

void DescriptiveIPTWithinHost::newInfection(){
  if (_MOI <= MAX_INFECTIONS) {
    _cumulativeInfections++;
    infections.push_back(new DescriptiveIPTInfection(_lastSPDose, Simulation::simulationTime));
    _MOI++;
  }
}

// -----    -----

void DescriptiveIPTWithinHost::clearInfections (bool isSevere) {
  int fortnight = int((14.0/Global::interval)+0.5);	// round to nearest
  if (isSevere) {
  } else if(Simulation::simulationTime-_lastIptiOrPlacebo <= fortnight) {
          // IPTi trials used quinine for fevers within 14 days of an ipti or placebo dose   
  } else if(Simulation::simulationTime-_lastSPDose <=  fortnight) {
          /*
    second line used if fever within 14 days of SP dose (ipti or treatment)
    TODO: if this code is to survive, then the iptiEffect values should be 
    symbolic constants
          */
  } else if( iptiEffect ==  2 ||  iptiEffect ==  12) {
    _lastSPDose=Simulation::simulationTime+1;
  } else if( iptiEffect ==  3 ||  iptiEffect ==  13) {
  } else if(iptiEffect >=  14 && iptiEffect < 30) {
  } else {
    _lastSPDose=Simulation::simulationTime+1;
    // SPAction will first act at the beginning of the next Global::interval
  }
  clearAllInfections();
}

void DescriptiveIPTWithinHost::IPTSetLastSPDose (int agetstep, int ageGroup) {
  if (Simulation::timeStep < 0) return;
  // assumes 5-day intervals and Niakhar seasonality
  // These numbers, should have MAX = MIN + 18 (modulo 73).
  static int IPT_MIN_INTERVAL[9] = { 43, 49, 55, 61, 67, 37, 31, 25, 19 };
  static int IPT_MAX_INTERVAL[9] = { 61, 67,  0,  6, 12, 55, 49, 43, 31 };
  
  if (iptiEffect >= 14 && iptiEffect <= 22) {
    int yearInterval = Simulation::simulationTime % Global::intervalsPerYear;
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
      if (gsl::rngUniform() <  iptiCoverage[i]) {
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

void DescriptiveIPTWithinHost::IPTiTreatment (int ageGroup) {
  //Set the last SP Dose given for the eligible humans - is this all we need to do?
  
  _lastIptiOrPlacebo = Simulation::simulationTime;
  
  // iptiEffect denotes treatment or placebo group
  // and also the treatment given when sick (trial-dependent)
  if (iptiEffect >= 10){
    _lastSPDose = Simulation::simulationTime;
    Simulation::gMainSummary->reportIPTDose(ageGroup);
  }
}


// -----  density calculation  -----

void DescriptiveIPTWithinHost::SPAction(){
  /*TODO if we want to look at presumptive SP treatment with the PkPD model we
  need to add some code here that will be conditionally implemented depending on the
  model version.*/

  std::list<DescriptiveInfection*>::iterator iter=infections.begin();
  while(iter != infections.end()) {
    if (1 + Simulation::simulationTime - (*iter)->getStartDate() > DescriptiveInfection::latentp) {
      DescriptiveIPTInfection* infec = dynamic_cast<DescriptiveIPTInfection*> (*iter);
      if (infec == NULL) throw logic_error ("infections should be of type DescriptiveInfection");
      //FIXME: put this if logic in DescriptiveIPTInfection and remove getGenoTypeID
      size_t genoTypeIndex = infec->getGenoTypeID()-1;
      if ((gsl::rngUniform() <= DescriptiveIPTInfection::genotypeACR[genoTypeIndex]) &&
	  (Simulation::simulationTime - _lastSPDose <= DescriptiveIPTInfection::genotypeProph[genoTypeIndex])) {
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

void DescriptiveIPTWithinHost::IPTattenuateAsexualMinTotalDensity () {
  //NOTE: the _cumulativeInfections>0 check is probably unintended, but was extracted from other logic and put here to preserve results.
  if (Global::modelVersion & ATTENUATION_ASEXUAL_DENSITY && _cumulativeInfections > 0) {
    if (_SPattenuationt > Simulation::simulationTime && totalDensity < 10) {
      totalDensity = 10;
      _cumulativeY += 10;
    }
  }
}

void DescriptiveIPTWithinHost::IPTattenuateAsexualDensity (DescriptiveInfection* inf) {
  if (!(Global::modelVersion & ATTENUATION_ASEXUAL_DENSITY)) return;
  
  DescriptiveIPTInfection* iptInf = dynamic_cast<DescriptiveIPTInfection*> (inf);
  if (iptInf == NULL)
    throw invalid_argument ("inf should be a DescriptiveIPTInfection");
  if (iptInf->doSPAttenuation()) {
    double attFact = iptInf->asexualAttenuation();
    timeStepMaxDensity *= attFact;
    _SPattenuationt = (int) std::max(double(_SPattenuationt),
				     iptInf->getAsexualAttenuationEndDate());
  }
}
