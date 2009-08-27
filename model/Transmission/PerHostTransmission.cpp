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
#include "summary.h"
#include "intervention.h"
#include "inputData.h"


// -----  PerHostTransmission static  -----

const double PerHostTransmission::bsa_prop[WithinHostModel::nages] = { 0.1843, 0.2225, 0.252, 0.2706, 0.2873, 0.3068, 0.3215, 0.3389, 0.3527, 0.3677, 0.3866, 0.3987, 0.4126, 0.4235, 0.441, 0.4564, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
double PerHostTransmission::ageSpecificRelativeAvailability[WithinHostModel::nages];
vector<double> PerHostTransmission::cntItnTargetAgeTStep;
vector<double> PerHostTransmission::cntItnCoverage;


void PerHostTransmission::initParameters (const scnXml::Interventions& interv) {
  for (size_t i=0; i<WithinHostModel::nages; i++) {
    ageSpecificRelativeAvailability[i] = bsa_prop[i] / (1-bsa_prop[i]);
  }
  
  if (interv.getContinuous().present()) {
    const scnXml::Continuous::ITNSequence& seqItn = interv.getContinuous().get().getITN();
    int n = seqItn.size();
    cntItnTargetAgeTStep.resize(n);
    cntItnCoverage.resize (n);
    
    for (int i=0;i<n; i++) {
      cntItnTargetAgeTStep[i] = static_cast<int>(floor(seqItn[i].getTargetAgeYrs() * daysInYear / (1.0*Global::interval)));
      cntItnCoverage[i] = seqItn[i].getCoverage();
    }
  }
}


// -----  PerHostTransmission non-static -----

PerHostTransmission::PerHostTransmission () :
    timestepITN(TIMESTEP_NEVER), timestepIRS(TIMESTEP_NEVER), timestepVA(TIMESTEP_NEVER),
    nextItnDistribution(0)
{}
void PerHostTransmission::initialise (TransmissionModel& tm, double availabilityFactor) {
  _entoAvailability = availabilityFactor;
  VectorTransmission* vTM = dynamic_cast<VectorTransmission*> (&tm);
  if (vTM) {
    species.resize (vTM->numSpecies);
    for (size_t i = 0; i < vTM->numSpecies; ++i)
      species[i].initialise (vTM->species[i].getHumanBase(), availabilityFactor);
  }
}

PerHostTransmission::PerHostTransmission (istream& in, TransmissionModel& tm) {
  in >> _entoAvailability;
  in >> timestepITN;
  in >> timestepIRS;
  in >> timestepVA;
  in >> nextItnDistribution;
  VectorTransmission* vTM = dynamic_cast<VectorTransmission*> (&tm);
  if (vTM) {
    species.resize (vTM->numSpecies);
    for (vector<HostMosquitoInteraction>::iterator hMI = species.begin(); hMI != species.end(); ++hMI)
      hMI->read (in);
  }
}

void PerHostTransmission::write (ostream& out) const {
  out << _entoAvailability << endl;
  out << timestepITN << endl;
  out << timestepIRS << endl;
  out << timestepVA << endl;
  out << nextItnDistribution << endl;
  for (vector<HostMosquitoInteraction>::const_iterator hMI = species.begin(); hMI != species.end(); ++hMI)
    hMI->write (out);
}

// NOTE: in the case an ITN / IRS is not present, this is only an approximation.
// But (Simulation::simulationTime - TIMESTEP_NEVER) is easily large enough for
// conceivable Weibull params that the value is 0.0 when rounded to a double.
// Performance-wise, an if() might give a small performance gain when
// interventions aren't present.
double PerHostTransmission::entoAvailabilityPartial (HostCategoryAnopheles& base, size_t speciesIndex) const {
  return species[speciesIndex].entoAvailability
    * (1.0 - base.ITNDeterrency (Simulation::simulationTime - timestepITN))
    * (1.0 - base.IRSDeterrency (Simulation::simulationTime - timestepIRS))
    * (1.0 - base.VADeterrency  (Simulation::simulationTime - timestepVA));
}
double PerHostTransmission::probMosqBiting (HostCategoryAnopheles& base, size_t speciesIndex) const {
  return species[speciesIndex].probMosqBiting
    * (1.0 - base.ITNPreprandialKillingEffect (Simulation::simulationTime - timestepITN));
}
double PerHostTransmission::probMosqResting (HostCategoryAnopheles& base, size_t speciesIndex) const {
  double P_C_i = species[speciesIndex].probMosqFindRestSite
    * (1.0 - base.ITNPostprandialKillingEffect (Simulation::simulationTime - timestepITN));
  double P_D_i = species[speciesIndex].probMosqSurvivalResting
    * (1.0 - base.IRSKillingEffect (Simulation::simulationTime - timestepIRS));
  return P_C_i * P_D_i;
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

void HostMosquitoInteraction::read (istream& in) {
  in >> entoAvailability;
  in >> probMosqBiting;
  in >> probMosqFindRestSite;
  in >> probMosqSurvivalResting;
}

void HostMosquitoInteraction::write (ostream& out) const {
  out << entoAvailability << endl;
  out << probMosqBiting << endl;
  out << probMosqFindRestSite << endl;
  out << probMosqSurvivalResting << endl;
}
