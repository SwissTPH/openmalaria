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

#ifndef Hmod_DrugModel
#define Hmod_DrugModel

#include <fstream>
using namespace std;

class ProteomeInstance;

/** Encapsulates both the static operations for drug models and the per-human
 * drug proxies.
 * 
 * Note that there currently needn't be a drug model, in which case an instance
 * of this class is created (which while inefficient, allows nicer code). */
class DrugModel {
public:
  ///@brief Static functions
  //@{
  static void init ();
  
  /// Read static variables from checkpoint
  static void readStatic (istream& in);
  /// Write static variables to checkpoint
  static void writeStatic (ostream& out);
  
  /// Create a new DrugModel
  static DrugModel* createDrugModel ();
  /// Load a DrugModel from a checkpoint
  static DrugModel* createDrugModel (istream& in);
  //@}
  
  /// @brief Constructors, destructor and checkpointing
  //@{
  /// Create a new instance
  DrugModel () {}
  /// Load an instance from a checkpoint
  DrugModel (istream& in) {}
  /// Destroy an instance
  virtual ~DrugModel () {}
  /// Write checkpoint
  virtual void write (ostream& out) {}
  //@}
  
  /** Medicate drugs to an individual, which act on infections the following
   * timesteps, until rendered ineffective by decayDrugs().
   *
   * \param drugAbbrev - The drug abbreviation.
   * \param qty        - the quantity (which units?).
   * \param time       - Time in minutes (or hours?) since start of the simulation tStep.
   *  \param weight    - Weight of human in kg */
  virtual void medicate(string drugAbbrev, double qty, int time, double weight) {}
  
  /// Called each timestep immediately after the drug acts on any infections.
  //NOTE: does calling after applying drug effects make the most sense for all models?
  virtual void decayDrugs () {}
  
  /** This is how drugs act on infections.
   *
   * Each timestep, on each infection, the parasite density is multiplied by
   * the return value of this infection. The WithinHostModels are responsible
   * for clearing infections once the parasite density is negligible. */
  virtual double getDrugFactor (ProteomeInstance* infProteome) {
    return 0.0;
  }
};

#endif
