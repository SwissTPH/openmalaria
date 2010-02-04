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

#ifndef Hmod_HoshenDrug
#define Hmod_HoshenDrug

#include <string>
#include <deque>
#include <map>
#include <vector>
#include "Dose.h"
#include "HoshenDrugType.h"
#include "Global.h"
#include "PkPd/Proteome.h"

using namespace std;

namespace OM { namespace PkPd {
    
/** A class holding hoshen pkpd drug use info.
 *
 * Each human has an instance for each type of drug present in their blood. */
class HoshenDrug {
public:
  /** Create a new instance. */
  HoshenDrug (const HoshenDrugType*);
  
  /** Called per timestep to reduce concentrations.
   *
   * If remaining concentration is negligible, return true, and this class
   * object will be deleted. */
  bool decay();
  
  /** Add a dose.
   *
   * @param concentration
   * @param time Days */
  void addDose (double concentration, double time);
  
  string getAbbreviation() const { return typeData->abbreviation;}
  //double getAbsorptionFactor() const { return typeData->absorptionFactor;}
  //double getHalfLife() const { return typeData->halfLife;}
  
  double getAbsorptionFactor() const { return ((HoshenDrugType*)typeData)->absorptionFactor;}

  double calculateDrugFactor(uint32_t proteome_ID) const;
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      _concentration & stream;
      _nextConcentration & stream;
  }
  
protected:
   /** Calculate multiplier to decay a concentration by a duration of time
    *
    * @param time Duration in days to decay over */
  double decayFactor (double time);
  
  /// Always links a drug instance to its drug-type data
  const HoshenDrugType* typeData;
  
  //! Drug concentration (ng/mL ?).
  double _concentration;
  //! Drug concentration on the next cycle (always should be whatever calcNextConcentration sets).
  double _nextConcentration;
};

} }
#endif