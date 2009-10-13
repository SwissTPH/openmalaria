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

#include "Global.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/Vector/VectorAnopheles.h"

namespace scnXml {
  class Vector;
}

//! Transmission models, Chitnis et al
class VectorTransmission : public TransmissionModel {
public:
  VectorTransmission(const scnXml::Vector vectorData);
  virtual ~VectorTransmission();
  
  virtual void writeV(ostream& out) const;
  virtual void readV(istream& in);
  
  /** Extra initialisation, requiring information from the human population
   * structure. */
  virtual void setupNv0 (const std::list<Human>& population, int populationSize);
  
  /** Length with which to force vector calculations, while waiting for human
   * population to stabilise. */
  virtual int vectorInitDuration () {
    // 30 years allows human availability & infectiousness to stabilise somewhat, at least according to one scenario.
    return Global::intervalsPerYear*30;
  }
  /** Called after end of vectorInitDuration() and after init iterations.
   *
   * Should determine whether another init iteration is needed, make necessary
   * adjustments, and return number of timesteps to run this initialisation for
   * (0 if a further iteration is not needed). */
  virtual int vectorInitIterate ();
  
  /** Initialise the main simulation. */
  void initMainSimulation ();
  
  virtual void vectorUpdate (const std::list<Human>& population, int simulationTime);
  
  /** Calculates EIR (in adults).
   * 
   * \param simulationTime Time since start of simulation . */
  virtual double calculateEIR(int simulationTime, PerHostTransmission& host, double ageInYears); 

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
  
  /** Per anopheles species data.
   *
   * Array will be recreated by constructor, but some members of VectorAnopheles
   * need to be checkpointed. */
  vector<VectorAnopheles> species;
  
  /** A map of anopheles species/variant name to an index in species.
   *
   * When the main ento data is read from XML, each anopheles section is read
   * into one index of the species array. The name and this index are added
   * here.
   * 
   * Other data read from XML should look up the name here and use the index
   * found. Doesn't need checkpointing. */
  map<string,size_t> speciesIndex;
  //@}
  
  friend class PerHostTransmission;
  friend class VectorAnophelesSuite;
};
#endif
