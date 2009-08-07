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

#include "Transmission/Vector.h"
#include "inputData.h"
#include <fstream>


/*****************************************************************************
 *                       New code, written in C++                            *
 *****************************************************************************/


VectorTransmission::VectorTransmission (const scnXml::Vector vectorData, const std::list<Human>& population, int populationSize) {
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
  
  species[0].initialise (anophelesList[0], 0, population, populationSize, initialisationEIR);
  simulationMode = species[0].getSimulationMode();
  for (size_t i = 1; i < numSpecies; ++i) {
    species[i].initialise (anophelesList[i], i, population, populationSize, initialisationEIR);
    
    if (simulationMode != species[i].getSimulationMode())
      // Probably fourier eir data present only for some species
      throw xml_scenario_error("Vector model requests EIR-driven initialisation for some species but not others");
  }
  
  // We want the EIR to effectively be the sum of the EIR for each day in the interval
  for (size_t i = 0; i < initialisationEIR.size(); ++i) {
    initialisationEIR[i] *= Global::interval;
  
    // Calculate total annual EIR
    annualEIR += initialisationEIR[i];
  }
}
VectorTransmission::~VectorTransmission () {
  for (size_t i = 0; i < numSpecies; ++i)
    species[i].destroy();
}

void VectorTransmission::initMainSimulation(const std::list<Human>& population, int populationSize) {
  cerr << "Warning: using incomplete VectorTransmission transmission model!" << endl;
  for (size_t i = 0; i < numSpecies; ++i)
    species[i].initMainSimulation (i, population, populationSize, kappa);
  simulationMode = dynamicEIR;	// currently always the case post-initialisation
}

double VectorTransmission::calculateEIR(int simulationTime, PerHostTransmission& host, double ageInYears) {
  double EIR = 0.0;
  for (size_t i = 0; i < numSpecies; ++i) {
    EIR += species[i].calculateEIR (i, host);
  }
  return EIR * PerHostTransmission::getRelativeAvailability(ageInYears);
}


// Every Global::interval days:
void VectorTransmission::advancePeriod (const std::list<Human>& population, int simulationTime) {
  if (simulationMode == equilibriumMode) return;
  for (size_t i = 0; i < numSpecies; ++i)
    species[i].advancePeriod (population, simulationTime, i);
}
