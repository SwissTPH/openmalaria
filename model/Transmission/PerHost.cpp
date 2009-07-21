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
#include "Transmission/PerHost.h"
#include "Transmission/Vector.h"
#include "summary.h"
#include "intervention.h"


// -----  PerHostTransmission static  -----

const double PerHostTransmission::bsa_prop[WithinHostModel::nages] = { 0.1843, 0.2225, 0.252, 0.2706, 0.2873, 0.3068, 0.3215, 0.3389, 0.3527, 0.3677, 0.3866, 0.3987, 0.4126, 0.4235, 0.441, 0.4564, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
double PerHostTransmission::ageSpecificRelativeAvailability[WithinHostModel::nages];


void PerHostTransmission::initParameters () {
  EntoInterventionITN::initParameters();
  EntoInterventionIRS::initParameters();
  
  for (size_t i=0; i<WithinHostModel::nages; i++) {
    ageSpecificRelativeAvailability[i] = bsa_prop[i] / (1-bsa_prop[i]);
  }
}


// -----  PerHostTransmission non-static -----

PerHostTransmission::PerHostTransmission () {
}
void PerHostTransmission::initialise (TransmissionModel& tm, double availabilityFactor) {
  _entoAvailability = availabilityFactor;
  VectorTransmission* vTM = dynamic_cast<VectorTransmission*> (&tm);
  if (vTM) {
    species.resize (vTM->numSpecies);
    for (size_t i = 0; i < vTM->numSpecies; ++i)
      species[i].initialise (vTM->species[i], availabilityFactor);
  }
}

PerHostTransmission::PerHostTransmission (istream& in, TransmissionModel& tm) {
  in >> _entoAvailability;
  VectorTransmission* vTM = dynamic_cast<VectorTransmission*> (&tm);
  if (vTM) {
    species.resize (vTM->numSpecies);
    for (vector<HostMosquitoInteraction>::iterator hMI = species.begin(); hMI != species.end(); ++hMI)
      hMI->read (in);
  }
}

void PerHostTransmission::write (ostream& out) const {
  out << _entoAvailability << endl;
  for (vector<HostMosquitoInteraction>::const_iterator hMI = species.begin(); hMI != species.end(); ++hMI)
    hMI->write (out);
}

//TODO: intervention effects on these parameters:
double PerHostTransmission::entoAvailabilityPartial (size_t speciesIndex) const {
  return species[speciesIndex].entoAvailability;
//   * species[speciesIndex].entoInterventionITN.availability()
//   * species[speciesIndex].entoInterventionIRS.availability();
}
double PerHostTransmission::probMosqBiting (size_t speciesIndex) const {
  return species[speciesIndex].probMosqBiting;
//   * species[speciesIndex].entoInterventionITN.probMosqBiting();
}
double PerHostTransmission::probMosqFindRestSite (size_t speciesIndex) const {
  return species[speciesIndex].probMosqFindRestSite;
//   * species[speciesIndex].entoInterventionITN.probMosqFindRestSite();
}
double PerHostTransmission::probMosqSurvivalResting (size_t speciesIndex) const {
  return species[speciesIndex].probMosqSurvivalResting;
//   * species[speciesIndex].entoInterventionIRS.probMosqSurvivalResting();
}


// ----- HostMosquitoInteraction non-static -----

void HostMosquitoInteraction::initialise (VectorTransmissionSpecies base, double availabilityFactor)
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
  in >> entoInterventionITN;
  in >> entoInterventionIRS;
}

void HostMosquitoInteraction::write (ostream& out) const {
  out << entoAvailability << endl;
  out << probMosqBiting << endl;
  out << probMosqFindRestSite << endl;
  out << probMosqSurvivalResting << endl;
  out << entoInterventionITN;
  out << entoInterventionIRS;
}
