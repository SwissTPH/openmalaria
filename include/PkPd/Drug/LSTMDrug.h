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

#ifndef Hmod_LSTMDrug
#define Hmod_LSTMDrug

#include "Drug.h"
#include "Dose.h"
#include "Global.h"
#include "PkPd/Proteome.h"
#include "LSTMDrugType.h"

using namespace std;

namespace OM { namespace PkPd {
    
/** A class holding pkpd drug use info.
 *
 * Each human has an instance for each type of drug present in their blood. */
class LSTMDrug : public Drug {
public:
  /** Create a new instance. */
  LSTMDrug (const LSTMDrugType*);
  
  string getAbbreviation() const { return typeData->abbreviation;}
  //double getAbsorptionFactor() const { return typeData->absorptionFactor;}
  //double getHalfLife() const { return typeData->halfLife;}
  
  /** Add amount to the concentration of drug, at time delay past the start of
   * the current timestep. */
  void storeDose (double amount, int delay);
  
  double calculateDrugFactor(const ProteomeInstance* infProteome, double ageYears, double weight_kg);
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      concentration & stream;
      doses & stream;
  }
  
protected:
  /** Calculate multiplier to decay a concentration by a duration of time
    *
    * @param time Duration in minutes to decay over */
  double decayFactor (double time);

  /// Always links a drug instance to its drug-type data
  const LSTMDrugType* typeData;
  
  double concentration;
  
  /// Per-dose information. Still to be properly defined.
  list<Dose> doses;
};

} }
#endif