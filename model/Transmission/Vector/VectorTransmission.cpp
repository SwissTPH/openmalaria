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

#include "Transmission/Vector/VectorTransmission.h"
#include "inputData.h"
#include "util/vectors.h"
#include "util/ModelOptions.hpp"

#include <fstream>

namespace OM { namespace Transmission {
    using namespace OM::util;

VectorTransmission::VectorTransmission (const scnXml::Vector vectorData)
{  
  for (size_t j=0;j<Global::intervalsPerYear; j++)
    initialisationEIR[j]=0.0;
  
  // Each item in the AnophelesSequence represents an anopheles species.
  // TransmissionModel::createTransmissionModel checks length of list >= 1.
  const scnXml::Vector::AnophelesSequence anophelesList = vectorData.getAnopheles();
  numSpecies = anophelesList.size();
  if (numSpecies < 1)
    throw util::xml_scenario_error ("Can't use Vector model without data for at least one anopheles species!");
#ifdef OMV_CSV_REPORTING
  csvReporting << "##\t##" << endl;	// live-graph needs a deliminator specifier when it's not a comma
  species.resize (numSpecies, VectorAnopheles(this,csvReporting));
  csvReporting << "simulation time\t";
#else
  species.resize (numSpecies, VectorAnopheles(this));
#endif
  
  for (size_t i = 0; i < numSpecies; ++i) {
    string name = species[i].initialise (anophelesList[i], i,
					 initialisationEIR);
    speciesIndex[name] = i;
    
#ifdef OMV_CSV_REPORTING
    csvReporting << "N_v0("<<i<<")\tN_v("<<i<<")\tO_v("<<i<<")\tS_v("<<i<<")\t";
#endif
  }
  
#ifdef OMV_CSV_REPORTING
  csvReporting << "input EIR\tresultant EIR\thuman infectiousness\thuman availability" << endl;
#endif
  
  for (size_t i = 0; i < initialisationEIR.size(); ++i) {
    // Calculate total annual EIR
    annualEIR += initialisationEIR[i];
  }
  
  
  // -----  Initialise interventions  -----
  const scnXml::Interventions::AnophelesSequence& intervSeq = InputData.getInterventions().getAnopheles();
  for (scnXml::Interventions::AnophelesSequence::const_iterator it = intervSeq.begin(); it != intervSeq.end(); ++it) {
    species[getSpeciesIndex(it->getMosquito())].setInterventionDescription (*it);
  }
  for (map<string,size_t>::const_iterator it = speciesIndex.begin(); it != speciesIndex.end(); ++it)
      species[it->second].checkInterventionDescriptions (it->first);
}
VectorTransmission::~VectorTransmission () {
  for (size_t i = 0; i < numSpecies; ++i)
    species[i].destroy();
}

void VectorTransmission::setupNv0 (const std::list<Host::Human>& population, int populationSize) {
  for (size_t i = 0; i < numSpecies; ++i) {
    species[i].setupNv0 (i, population, populationSize);
  }
}

int VectorTransmission::vectorInitIterate () {
  bool iterate = false;
  for (size_t i = 0; i < numSpecies; ++i)
    iterate |= species[i].vectorInitIterate ();
  if (iterate) {
    simulationMode = equilibriumMode;
    return Global::intervalsPerYear*2;	//TODO: how long?
  } else {
    simulationMode = dynamicEIR;
    return 0;
  }
}

void VectorTransmission::initMainSimulation() {
  // Check every time at end of init that, to a very low tolerence, the average
  // EIR produced is what was expected. Note: it's become clear the simulated
  // EIR is never going to very accurately match the input EIR.
  if (!vectors::approxEqual(initialisationEIR, innoculationsPerDayOfYear, 1)) {
    cerr << "Generated EIR not as expected (expected, generated):\n";
    cerr << initialisationEIR << '\n' << innoculationsPerDayOfYear << endl;
  }
  
  simulationMode = InputData.get_mode();	// allow forcing equilibrium mode like with non-vector model
  if (simulationMode != 2 && simulationMode != 4)
    throw util::xml_scenario_error("mode attribute has invalid value (expected: 2 or 4)");
}

double VectorTransmission::calculateEIR(int simulationTime, PerHostTransmission& host, double ageInYears) {
  if (simulationMode == equilibriumMode)
    return initialisationEIR[simulationTime%Global::intervalsPerYear]
	 * host.relativeAvailabilityHetAge (ageInYears) * ageCorrectionFactor;
  
  double EIR = 0.0;
  for (size_t i = 0; i < numSpecies; ++i) {
    EIR += species[i].calculateEIR (i, host);
  }
  return EIR * PerHostTransmission::relativeAvailabilityAge (ageInYears) * ageCorrectionFactor;
}


// Every Global::interval days:
void VectorTransmission::vectorUpdate (const std::list<Host::Human>& population, int simulationTime) {
#ifdef OMV_CSV_REPORTING
  csvReporting << simulationTime << '\t';
#endif
  for (size_t i = 0; i < numSpecies; ++i)
    species[i].advancePeriod (population, simulationTime, i, simulationMode == dynamicEIR);
}

void VectorTransmission::intervLarviciding (const scnXml::Larviciding& anoph) {
  const scnXml::Larviciding::AnophelesSequence& seq = anoph.getAnopheles();
  for (scnXml::Larviciding::AnophelesSequence::const_iterator it = seq.begin(); it != seq.end(); ++it)
    species[getSpeciesIndex(it->getMosquito())].intervLarviciding(*it);
}

void VectorTransmission::summarize (Survey& survey) {
    TransmissionModel::summarize (survey);
    
    for (map<string,size_t>::const_iterator it = speciesIndex.begin(); it != speciesIndex.end(); ++it)
	species[it->second].summarize (it->first, survey);
}


void VectorTransmission::checkpoint (istream& stream) {
    TransmissionModel::checkpoint (stream);
#ifdef OMV_CSV_REPORTING
    util::checkpoint::checkpoint (species, stream, VectorAnopheles (this, csvReporting));
#else
    util::checkpoint::checkpoint (species, stream, VectorAnopheles (this));
#endif
}
void VectorTransmission::checkpoint (ostream& stream) {
    TransmissionModel::checkpoint (stream);
    species & stream;
}

} }