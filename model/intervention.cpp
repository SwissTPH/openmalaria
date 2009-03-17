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

#include "intervention.h"
#include "inputData.h"
#include "GSLWrapper.h"
#include <iostream>

enum VaccineType {
  preerythrocytic_reduces_h = 1,
  erythrocytic_reduces_y = 2,
  transmission_blocking_reduces_k= 3,
};

bool Vaccine::anyVaccine = false;
double *Vaccine::targetagetstep;
double *Vaccine::vaccineCoverage;
int Vaccine::_numberOfEpiDoses = 0;
Vaccine Vaccine::PEV;
Vaccine Vaccine::BSV;
Vaccine Vaccine::TBV;

short IPTIntervention::IPT;
int IPTIntervention::numberOfIPTiDoses;
int *IPTIntervention::iptiTargetagetstep;
double *IPTIntervention::iptiCoverage;
int IPTIntervention::iptiEffect;
int IPTIntervention::numberOfGenoTypes;
double *IPTIntervention::genotypeFreq;
double *IPTIntervention::genotypeACR;
int *IPTIntervention::genotypeProph;
int *IPTIntervention::genotypeTolPeriod;
double *IPTIntervention::genotypeAtten;


double Vaccine::getEfficacy (int numPrevDoses) {
  /* If initialMeanEfficacy.size or more doses have already been given, use
   * the last efficacy. */
  if (numPrevDoses >= initialMeanEfficacy.size())
    numPrevDoses = initialMeanEfficacy.size() - 1;
  if (initialMeanEfficacy[numPrevDoses] <  1) {
    double a = efficacyB * initialMeanEfficacy[numPrevDoses] / (1.0-initialMeanEfficacy[numPrevDoses]);
    return W_BETA(a, efficacyB);
  }
  else
    return 1.0;
}

void Vaccine::initParameters() {
  const VaccineDescription *VdPEV = 0, *VdBSV = 0, *VdTBV = 0;
  const Interventions& interventions = getInterventions();
  const Interventions::VaccineDescriptionSequence& vaccDesc = interventions.getVaccineDescription();
  for (Interventions::VaccineDescriptionConstIterator i = vaccDesc.begin();
       i != vaccDesc.end(); i++) {
    int type = i->getVaccineType();
    if (type == preerythrocytic_reduces_h)
      VdPEV = &(*i);
    else if (type == erythrocytic_reduces_y)
      VdBSV = &(*i);
    else if (type == transmission_blocking_reduces_k)
      VdTBV = &(*i);
    
    anyVaccine = true;
  }
  if (!anyVaccine) return;

  //Read in vaccine specifications
  PEV.initVaccine (VdPEV);
  BSV.initVaccine (VdBSV);
  TBV.initVaccine (VdTBV);
  
  if (interventions.getContinuous().present())
    _numberOfEpiDoses = interventions.getContinuous().get().getVaccine().size();
  if (_numberOfEpiDoses) {
    targetagetstep = (double*)malloc(((_numberOfEpiDoses))*sizeof(double));
    vaccineCoverage = (double*)malloc(((_numberOfEpiDoses))*sizeof(double));
    const Continuous::VaccineSequence& cVS = interventions.getContinuous().get().getVaccine();
    for (int i=0;i<_numberOfEpiDoses; i++) {
      if (i >= cVS.size()) {
        cerr << "Expected " << _numberOfEpiDoses << " vaccine parameters in scenario.xml: interventions->continuous" << endl;
        throw 0;
      }
      targetagetstep[i] = floor(cVS[i].getTargetAgeYrs() * daysInYear/(double)Global::interval);
      vaccineCoverage[i] = cVS[i].getCoverage();
    }
  }
}

void Vaccine::initVaccine (const VaccineDescription* vd) {
  if (vd != NULL) {
    active = true;
    
    // set efficacyB:
    efficacyB = vd->getEfficacyB().getValue();
    
    // set initialMeanEfficacy:
    const VaccineDescription::InitialEfficacySequence ies = vd->getInitialEfficacy();
    initialMeanEfficacy.resize (ies.size());
    for (int i = 0; i < initialMeanEfficacy.size(); ++i)
      initialMeanEfficacy[i] = ies[i].getValue();
    
    // now, if halfLifeYrs > 0, calculate delay:
    double halfLifeYrs = vd->getHalfLifeYrs().getValue();
    if (halfLifeYrs <= 0)
      decay = 1.0;
    else
      decay = exp(-log(2.0) / (halfLifeYrs * daysInYear / (double)Global::interval));
  }
}

void Vaccine::clearParameters () {
  if (!Vaccine::anyVaccine)
    return;
  
  if (_numberOfEpiDoses == 0)
    return;
  free(targetagetstep);
  free(vaccineCoverage);
}

void IPTIntervention::initParameters () {
  const Interventions& xmlInterventions = getInterventions();
  IPT = xmlInterventions.getIptiDescription().present();
  if (!IPT)
    return;
  
  // --- IptiDescription begin ---
  const IptDescription& xmlIPTI = xmlInterventions.getIptiDescription().get();
  
  const IptDescription::InfGenotypeSequence& genotypes = xmlIPTI.getInfGenotype();
  numberOfGenoTypes = genotypes.size();
  genotypeFreq	= (double*)malloc(((numberOfGenoTypes))*sizeof(double));
  genotypeACR	= (double*)malloc(((numberOfGenoTypes))*sizeof(double));
  genotypeProph	= (int*)malloc(((numberOfGenoTypes))*sizeof(int));
  genotypeTolPeriod = (int*)malloc(((numberOfGenoTypes))*sizeof(int));
  genotypeAtten	= (double*)malloc(((numberOfGenoTypes))*sizeof(double));
  
  size_t i = 0;
  for (IptDescription::InfGenotypeConstIterator it = genotypes.begin(); it != genotypes.end(); ++it, ++i) {
    genotypeFreq[i]	= it->getFreq();
    genotypeACR[i]	= it->getACR();
    genotypeProph[i]	= it->getProph();
    genotypeTolPeriod[i]= it->getTolPeriod();
    genotypeAtten[i]	= it->getAtten();
  }
  
  iptiEffect = xmlIPTI.getIptiEffect();
  // --- IptiDescription end ---
  
  if (xmlInterventions.getContinuous().present()) {
    const Continuous::IptiSequence& xmlIpti = xmlInterventions.getContinuous().get().getIpti();
    numberOfIPTiDoses = xmlIpti.size();
    
    iptiTargetagetstep = (int*)malloc(((numberOfIPTiDoses))*sizeof(int));
    iptiCoverage = (double*)malloc(((numberOfIPTiDoses))*sizeof(double));
    for (int i=0;i<numberOfIPTiDoses; i++) {
      iptiTargetagetstep[i] = floor(xmlIpti[i].getTargetAgeYrs() * daysInYear / (1.0*Global::interval));
      iptiCoverage[i] = xmlIpti[i].getCoverage();
    }
  } else
    numberOfIPTiDoses = 0;
}

void IPTIntervention::clearParameters () {
  if (!IPT)
    return;
  
  free(genotypeFreq);
  free(genotypeACR);
  free(genotypeProph);
  free(genotypeTolPeriod);
  free(genotypeAtten);
  if (numberOfIPTiDoses) {
    free(iptiTargetagetstep);
    free(iptiCoverage);
  }
}
