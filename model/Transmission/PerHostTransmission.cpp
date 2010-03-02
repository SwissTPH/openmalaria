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
#include "Transmission/PerHostTransmission.h"
#include "Transmission/Vector/VectorTransmission.h"	//TODO: remove dependency
#include "inputData.h"

namespace OM { namespace Transmission {
        using namespace OM::util;

// -----  PerHostTransmission static  -----

vector<double> PerHostTransmission::cntItnTargetAgeTStep;
vector<double> PerHostTransmission::cntItnCoverage;


void PerHostTransmission::initParameters (const scnXml::Interventions& interv) {
    AgeGroupData::initParameters ();
  if (interv.getContinuous().present()) {
    const scnXml::Continuous::ITNSequence& seqItn = interv.getContinuous().get().getITN();
    int n = seqItn.size();
    cntItnTargetAgeTStep.resize(n);
    cntItnCoverage.resize (n);
    
    for (int i=0;i<n; i++) {
      cntItnTargetAgeTStep[i] = static_cast<int>(floor(seqItn[i].getTargetAgeYrs() * Global::DAYS_IN_YEAR / (1.0*Global::interval)));
      cntItnCoverage[i] = seqItn[i].getCoverage();
    }
  }
}


// -----  PerHostTransmission non-static -----

PerHostTransmission::PerHostTransmission () :
    timestepITN(Global::TIMESTEP_NEVER), timestepIRS(Global::TIMESTEP_NEVER), timestepVA(Global::TIMESTEP_NEVER),
    nextItnDistribution(0)
{}
void PerHostTransmission::initialise (TransmissionModel& tm, double availabilityFactor) {
  _relativeAvailabilityHet = availabilityFactor;
  VectorTransmission* vTM = dynamic_cast<VectorTransmission*> (&tm);
  if (vTM) {
    species.resize (vTM->numSpecies);
    for (size_t i = 0; i < vTM->numSpecies; ++i)
      species[i].initialise (vTM->species[i].getHumanBase(), availabilityFactor);
  }
}


// Note: in the case an intervention is not present, we can use the approximation
// of Weibull decay over (Global::simulationTime - Global::TIMESTEP_NEVER) timesteps
// (easily large enough for conceivable Weibull params that the value is 0.0 when
// rounded to a double. Performance-wise it's perhaps slightly slower than using
// an if() when interventions aren't present.
double PerHostTransmission::entoAvailabilityHetVecItv (const HostCategoryAnopheles& base, size_t speciesIndex) const {
  double alpha_i = species[speciesIndex].entoAvailability;
  if (timestepITN >= 0)
    alpha_i *= (1.0 - base.ITNDeterrency (Global::simulationTime - timestepITN));
  if (timestepIRS >= 0)
    alpha_i *= (1.0 - base.IRSDeterrency (Global::simulationTime - timestepIRS));
  if (timestepVA >= 0)
    alpha_i *= (1.0 - base.VADeterrency  (Global::simulationTime - timestepVA));
  return alpha_i;
}
double PerHostTransmission::probMosqBiting (const HostCategoryAnopheles& base, size_t speciesIndex) const {
  double P_B_i = species[speciesIndex].probMosqBiting;
  if (timestepITN >= 0)
    P_B_i *= (1.0 - base.ITNPreprandialKillingEffect (Global::simulationTime - timestepITN));
  return P_B_i;
}
double PerHostTransmission::probMosqResting (const HostCategoryAnopheles& base, size_t speciesIndex) const {
  double P_C_i = species[speciesIndex].probMosqFindRestSite;
  if (timestepITN >= 0)
    P_C_i *= (1.0 - base.ITNPostprandialKillingEffect (Global::simulationTime - timestepITN));
  double P_D_i = species[speciesIndex].probMosqSurvivalResting;
  if (timestepIRS >= 0)
    P_D_i *= (1.0 - base.IRSKillingEffect (Global::simulationTime - timestepIRS));
  return P_C_i * P_D_i;
}

void PerHostTransmission::continousItnDistribution (int ageTSteps) {
  if (Global::timeStep >= 0 && nextItnDistribution < cntItnTargetAgeTStep.size()
    && cntItnTargetAgeTStep[nextItnDistribution] == ageTSteps) {
    if (random::uniform_01() < cntItnCoverage[nextItnDistribution])
      setupITN ();
    ++nextItnDistribution;
  }
}


// ----- HostMosquitoInteraction non-static -----

void HostMosquitoInteraction::initialise (HostCategoryAnopheles& base, double availabilityFactor)
{
  //TODO: vary to simulate heterogeneity
  
  entoAvailability = base.entoAvailability * availabilityFactor;
  probMosqBiting = base.probMosqBiting;
  probMosqFindRestSite = base.probMosqFindRestSite;
  probMosqSurvivalResting = base.probMosqSurvivalResting;
}

} }