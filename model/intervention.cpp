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

enum VaccineType {
  preerythrocytic_reduces_h = 1,
  erythrocytic_reduces_y = 2,
  transmission_blocking_reduces_k= 3,
};

bool Vaccine::anyVaccine = false;
double *Vaccine::targetagetstep;
double *Vaccine::vaccineCoverage;
size_t Vaccine::_numberOfEpiDoses = 0;
Vaccine Vaccine::PEV;
Vaccine Vaccine::BSV;
Vaccine Vaccine::TBV;


double Vaccine::getEfficacy (int numPrevDoses) {
  /* If initialMeanEfficacy.size or more doses have already been given, use
   * the last efficacy. */
  if (numPrevDoses >= (int)initialMeanEfficacy.size())
    numPrevDoses = initialMeanEfficacy.size() - 1;
  if (initialMeanEfficacy[numPrevDoses] <  1) {
    double a = efficacyB * initialMeanEfficacy[numPrevDoses] / (1.0-initialMeanEfficacy[numPrevDoses]);
    return W_BETA(a, efficacyB);
  }
  else
    return 1.0;
}

void Vaccine::initParameters() {
  const scnXml::VaccineDescription *VdPEV = 0, *VdBSV = 0, *VdTBV = 0;
  const scnXml::Interventions& interventions = getInterventions();
  const scnXml::Interventions::VaccineDescriptionSequence& vaccDesc = interventions.getVaccineDescription();
  for (scnXml::Interventions::VaccineDescriptionConstIterator i = vaccDesc.begin();
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
    const scnXml::Continuous::VaccineSequence& cVS = interventions.getContinuous().get().getVaccine();
    for (size_t i=0;i<_numberOfEpiDoses; i++) {
      if (i >= cVS.size()) {
	ostringstream msg;
	msg << "Expected " << _numberOfEpiDoses << " vaccine parameters in scenario.xml: interventions->continuous";
        throw xml_scenario_error (msg.str());
      }
      targetagetstep[i] = floor(cVS[i].getTargetAgeYrs() * daysInYear/(double)Global::interval);
      vaccineCoverage[i] = cVS[i].getCoverage();
    }
  }
}

void Vaccine::initVaccine (const scnXml::VaccineDescription* vd) {
  if (vd != NULL) {
    active = true;
    
    // set efficacyB:
    efficacyB = vd->getEfficacyB().getValue();
    
    // set initialMeanEfficacy:
    const scnXml::VaccineDescription::InitialEfficacySequence ies = vd->getInitialEfficacy();
    initialMeanEfficacy.resize (ies.size());
    for (size_t i = 0; i < initialMeanEfficacy.size(); ++i)
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
