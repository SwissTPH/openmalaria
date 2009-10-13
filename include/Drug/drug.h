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

#ifndef Hmod_drug
#define Hmod_drug

#include <string>
#include <deque>
#include <map>
#include <vector>
#include "Global.h"
#include "proteome.h"

using namespace std;


/** Information about each (type of) drug (rather than each use of a drug).
 *
 * Static data contains a list of all available drug types.
 * 
 * No DrugType data is checkpointed, because it is loaded by init() from XML
 * data. (Although if it cannot be reproduced by reloading it should be
 * checkpointed.) */
class DrugType {
public:
  ///@brief Static functions
  //@{
  /** Initialise the drug model. Called at start of simulation. */
  static void init ();
  
  //! Adds a new drug to the list
  static void addDrug(DrugType* drug);

  /** Find a DrugType by its abbreviation, and create a new Drug from that.
   *
   * Throws if the drug isn't found, so you can rely on it returning a valid
   * drug if it returns (doesn't throw). */
  static const DrugType* getDrug(string abbreviation);
  //@}
  
  ///@brief Non static (per instance) functions
  //@{
  /** Create a new DrugType.
   *
   * @param name	Name of the drug
   * @param abbreviation	Abbreviated name (e.g. CQ)
   * @param absorptionFactor	
   * @param halfLife	Half life of decay, in minutes
   */
  DrugType (string name, string abbreviation, double absorptionFactor, double halfLife);
  ~DrugType ();
  /* Checkpointing functions, which we shouldn't need now. If they are needed:
  /// Load an instance from a checkpoint.
  DrugType (istream& in);
  /// Write instance data to a checkpoint.
  void write (ostream& out) const;
  */

  //! Adds a PD Rule.
  /*! The order of rule adding is important! The first add should be the
   *  one with most mutations (typically the most resistant), the last
   *  one should be the sensitive (ie vector<mutation> = 0).
   */
  void addPDRule(vector<Mutation*> requiredMutations, double pdFactor);

  //! Parses the proteme instances.
  /*! Creates an association between ProteomeInstance and PD factor.
   *  This is solely for performance purposes.
   */
  void parseProteomeInstances();
  //@}
  
private:
  // The list of available drugs. Not checkpointed; should be set up by init().
  static map<string,DrugType> available;
  
  //BEGIN Drug-type fields (same for all drugs of same type)
  //! The drug abbreviated name, used for registry lookups.
  string abbreviation;
  //! The drug name.
  string name; 
  //! Absorption factor.
  /*! Absorption = dose * factor / weight
   */
  double absorptionFactor;
  //! Half-life (in minutes)
  double halfLife;
  //! Pharma dynamic list of parameters.
  /*! A ordered list of required mutations.
   * The parameter value can be found on pdParameters.
   * The order is important, the first one takes precedence
   *   (a map cannot impement this).
   */
  vector< vector<Mutation*> > requiredMutations;
  //! PD parameters (check requiredMutations)
  vector<double> pdParameters;
  //! Fast data structure to know the PD param per proteome
  map<int,double> proteomePDParameters;
  //END
  
  // Allow the Drug class to access private members
  friend class Drug;
};

/** A simple class to hold dose info. */
class Dose {
public:
  /** Create a new dose. */
  Dose (double x, double y) {
    this->x = x;
    this->y = y;
  }
  /** Load from a checkpoint. */
  Dose (istream& in);
  /** Write a checkpoint. */
  void write (ostream& out) const;
  
  /// Some type of data is wanted... (concentration at start of next timestep, and integral of concentration for this timestep?)
  double x,y;
};

/** A class holding drug use info.
 *
 * Each human has an instance for each type of drug present in their blood. */
class Drug {
public:
  /** Initialise the drug model. Called at start of simulation. */
  static void init ();
  
  /** Create a new instance. */
  Drug (const DrugType*);
  /** Load an instance from a checkpoint. */
  Drug (const DrugType*, istream& in);
  /** Write instance data to a checkpoint. */
  void write (ostream& out) const;
  
  string getAbbreviation() const { return typeData->abbreviation;}
  double getAbsorptionFactor() const { return typeData->absorptionFactor;}
  //double getHalfLife() const { return typeData->halfLife;}
  
  /** Add amount to the concentration of drug, at time delay past the start of
   * the current timestep. */
  void addDose (double amount, int delay);
  
  double getConcentration() const { return _concentration;}
  double getNextConcentration() const { return _nextConcentration;}
  double calculateDrugFactor(ProteomeInstance* infProteome) const;
  /** Called per timestep to reduce concentrations.
   *
   * If remaining concentration is negligible, return true, and this class
   * object will be deleted. */
  bool decay();
  
private:
  /** Calculate multiplier to decay a concentration by a duration of time
   *
   * @param time Duration in minutes to decay over */
  double decayFactor (double time);
  
  static double minutesPerTimeStep;
  
  /// Always links a drug instance to its drug-type data
  const DrugType* typeData;
  
  //! Drug concentration (ng/mL ?).
  double _concentration;
  //! Drug concentration on the next cycle (always should be whatever calcNextConcentration sets).
  double _nextConcentration;
  
  /// Per-dose information. Still to be properly defined.
  deque<Dose> doses;
};

#endif
