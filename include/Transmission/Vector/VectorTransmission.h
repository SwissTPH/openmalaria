/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef Hmod_VectorTransmission
#define Hmod_VectorTransmission
/* This file should contain the headers of all routines that we write in the C
 * program.
 */ 

/* We also include library headers here. */ 
#include "Transmission/TransmissionModel.h"
#include "Transmission/Vector/VectorAnopheles.h"

//#define VectorTransmission_PRINT_CalcInitMosqEmergeRate

namespace scnXml {
  class Vector;
}

//! Transmission models, Chitnis et al
class VectorTransmission : public TransmissionModel {
public:
  VectorTransmission(const scnXml::Vector vectorData);
  virtual ~VectorTransmission();
  
  /** Extra initialisation, requiring information from the human population
   * structure. */
  virtual void setupNv0 (const std::list<Human>& population, int populationSize);
  
  /** Length with which to force input EIR prior to switching to dynamic EIR
   * and starting updateOneLifespan(). Note: no explicit reason why these
   * events should coincide. */
  virtual int vectorInitDuration () {
    // NOTE: One year should be ample time to infect a significant number of
    // humans and start the infection cycle, but probably not long enough to
    // reach an equilibrium state.
    return Global::intervalsPerYear;
  }
  virtual void endVectorInitPeriod ();
  
  /** Initialise the main simulation.
   * 
   * Calculates mosquito emergence rate.
   * 
   * \param populationSize The total number of hosts.
   * 
   * Emergence rate calculations assume only one type of host; i.e. it calculates
   * the rate for a stable situation before interventions are introduced. */
  void initMainSimulation (const std::list<Human>& population, int populationSize); 

  /** Calculates EIR (in adults).
   * 
   * \param simulationTime Time since start of simulation . */
  virtual double calculateEIR(int simulationTime, PerHostTransmission& host, double ageInYears); 

  /** This needs to be called every interval. */
  virtual void advancePeriod (const std::list<Human>& population, int simulationTime);
  
  virtual void intervLarviciding (const scnXml::Larviciding&);
  
private:
  /** Return the index in speciesIndex of mosquito, throwing if not found. */
  size_t getSpeciesIndex (string mosquito) {
    map<string,size_t>::const_iterator sIndex = speciesIndex.find (mosquito);
    if (sIndex == speciesIndex.end()) {
      ostringstream oss;
      oss << "Intervention description for unincluded anopheles species \""
	  << mosquito << '"';
      throw xml_scenario_error(oss.str());
    }
    return sIndex->second;
  }
  
  /** @brief Access to per (anopheles) species data.
   *
   * Set by constructor so don't checkpoint. */
  //@{
  /** The number of discrete species of anopheles mosquitos to be modelled.
   *
   * Must be the same as species.size(). */
  size_t numSpecies;
  
  /** Per anopheles species data. */
  vector<VectorAnopheles> species;
  
  /** A map of anopheles species/variant name to an index in species.
   *
   * When the main ento data is read from XML, each anopheles section is read
   * into one index of the species array. The name and this index are added
   * here.
   * 
   * Other data read from XML should look up the name here and use the index
   * found. */
  map<string,size_t> speciesIndex;
  //@}
  
  /*NOTE: add NonHumanHosts data here:
  per-species parameters
  number of hosts */
  
  friend class PerHostTransmission;
  friend class VectorAnophelesSuite;
};
#endif
