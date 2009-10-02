/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

// Variable names largely come from Nakul Chitnis's paper:
// "A mathematical model for the dynamics of malaria in mosquitoes feeding on
// a heterogeneous host population" (3rd Oct. 2007).

/* Entomology model coordinator: Nakul Chitnis. */ 

#include "Transmission/Vector/VectorTransmission.h"
#include "inputData.h"
#include "util/vectors.h"
#include <fstream>


/*****************************************************************************
 *                       New code, written in C++                            *
 *****************************************************************************/


VectorTransmission::VectorTransmission (const scnXml::Vector vectorData) {
  if ((Global::modelVersion & (NEGATIVE_BINOMIAL_MASS_ACTION|LOGNORMAL_MASS_ACTION))==0)
    throw xml_scenario_error ("VectorTransmission is incompatible with the original InfectionIncidenceModel");
  
  for (size_t j=0;j<Global::intervalsPerYear; j++)
    initialisationEIR[j]=0.0;
  
  // Each item in the AnophelesSequence represents an anopheles species.
  // TransmissionModel::createTransmissionModel checks length of list >= 1.
  const scnXml::Vector::AnophelesSequence anophelesList = vectorData.getAnopheles();
  numSpecies = anophelesList.size();
  if (numSpecies < 1)
    throw xml_scenario_error ("Can't use Vector model without data for at least one anopheles species!");
  species.resize (numSpecies);
  
  for (size_t i = 0; i < numSpecies; ++i) {
    string name = species[i].initialise (anophelesList[i], i,
					 initialisationEIR);
    speciesIndex[name] = i;
  }
  
  // We want the EIR to effectively be the sum of the EIR for each day in the interval
  for (size_t i = 0; i < initialisationEIR.size(); ++i) {
    initialisationEIR[i] *= Global::interval;
  
    // Calculate total annual EIR
    annualEIR += initialisationEIR[i];
  }
  
  
  // -----  Initialise interventions  -----
  const scnXml::Interventions::AnophelesSequence& intervSeq = getInterventions().getAnopheles();
  for (scnXml::Interventions::AnophelesSequence::const_iterator it = intervSeq.begin(); it != intervSeq.end(); ++it) {
    species[getSpeciesIndex(it->getMosquito())].setInterventionDescription (*it);
  }
}
VectorTransmission::~VectorTransmission () {
  for (size_t i = 0; i < numSpecies; ++i)
    species[i].destroy();
}

void VectorTransmission::setupNv0 (const std::list<Human>& population, int populationSize) {
  for (size_t i = 0; i < numSpecies; ++i) {
    species[i].setupNv0 (i, population, populationSize);
  }
}

void VectorTransmission::endVectorInitPeriod () {
  simulationMode = dynamicEIR;
}

void VectorTransmission::initMainSimulation() {
  // Check every time at end of init that, to a low tolerence,
  // the average EIR produced is what was expected:
  if (!vectors::approxEqual(initialisationEIR, innoculationsPerDayOfYear)) {
    cerr << "Generated EIR not as expected (expected, generated):\n";
    cerr << initialisationEIR << '\n' << innoculationsPerDayOfYear << endl;
  }
  
  simulationMode = get_mode();	// allow forcing equilibrium mode like with non-vector model
  if (simulationMode != 2 && simulationMode != 4)
    throw xml_scenario_error("mode attribute has invalid value (expected: 2 or 4)");
}

double VectorTransmission::calculateEIR(int simulationTime, PerHostTransmission& host, double ageInYears) {
  double EIR = 0.0;
  for (size_t i = 0; i < numSpecies; ++i) {
    EIR += species[i].calculateEIR (i, host);
  }
  return EIR * PerHostTransmission::getRelativeAvailability(ageInYears);
}


// Every Global::interval days:
void VectorTransmission::advanceStepCalcs (const std::list<Human>& population, int simulationTime, double& sumWeight, double& sumWt_kappa) {
  for (size_t i = 0; i < numSpecies; ++i)
    species[i].advancePeriod (population, simulationTime, i, simulationMode == dynamicEIR, sumWeight, sumWt_kappa, kappaByAge, nByAge);
}

void VectorTransmission::intervLarviciding (const scnXml::Larviciding& anoph) {
  const scnXml::Larviciding::AnophelesSequence& seq = anoph.getAnopheles();
  for (scnXml::Larviciding::AnophelesSequence::const_iterator it = seq.begin(); it != seq.end(); ++it)
    species[getSpeciesIndex(it->getMosquito())].intervLarviciding(*it);
}


void VectorTransmission::writeV (ostream& out) const {
  for (size_t i = 0; i < numSpecies; ++i)
    species[i].write (out);
}
void VectorTransmission::readV (istream& in) {
  for (size_t i = 0; i < numSpecies; ++i)
    species[i].read (in);
}
