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

#include <string>
#include <deque>
#include <map>
#include <vector>
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
  /** Load an instance from a checkpoint. */
  LSTMDrug (const LSTMDrugType*, istream& in);
  void write (ostream& out) const;
  
  /** Add amount to the concentration of drug, at time delay past the start of
   * the current timestep. */
  void addDose (double amount, int delay);
  
  virtual double calculateDrugFactor(const ProteomeInstance* infProteome) const;
  double getAbsorptionFactor() const { return ((LSTMDrugType*)typeData)->absorptionFactor;}


protected:
  /** Calculate multiplier to decay a concentration by a duration of time
    *
    * @param time Duration in minutes to decay over */
  virtual double decayFactor (double time);

  
  /// Per-dose information. Still to be properly defined.
  deque<Dose> doses;
};

} }
#endif