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

#ifndef Hmod_Drug
#define Hmod_Drug

#include <string>
#include <deque>
#include <map>
#include <vector>
#include "DrugType.h"
#include "Global.h"
#include "proteome.h"

using namespace std;


/** A class holding drug use info. This is an abstract base class, so it doesn't
 * include all details required.
 *
 * Each human has an instance for each type of drug present in their blood. */
class Drug {
public:
  /** Initialise the drug model. Called at start of simulation. */
  static void init ();
  
  
  string getAbbreviation() const { return typeData->abbreviation;}
  double getAbsorptionFactor() const { return typeData->absorptionFactor;}
  //double getHalfLife() const { return typeData->halfLife;}
  
  /** Add amount to the concentration of drug at time the start of the current
   * timestep (delay is expected to be 0). */
  void addDose (double amount, int delay);
  
  virtual double calculateDrugFactor(const ProteomeInstance* infProteome) const =0;
  /** Called per timestep to reduce concentrations.
   *
   * If remaining concentration is negligible, return true, and this class
   * object will be deleted. */
  bool decay();
  
protected:
  /** Create a new instance. */
  Drug (const DrugType*);
  /** Load an instance from a checkpoint. */
  Drug (const DrugType*, istream& in);
  /** Write instance data to a checkpoint. */
  void write (ostream& out) const;
  /** Calculate multiplier to decay a concentration by a duration of time
   *
   * @param time Duration in minutes to decay over */
  
  virtual double decayFactor (double time) =0;
  
  static double minutesPerTimeStep;
  
  /// Always links a drug instance to its drug-type data
  const DrugType* typeData;
  
  //! Drug concentration (ng/mL ?).
  double _concentration;
  //! Drug concentration on the next cycle (always should be whatever calcNextConcentration sets).
  double _nextConcentration;
};

#endif
