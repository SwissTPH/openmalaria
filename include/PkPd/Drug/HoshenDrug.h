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
#include "Drug.h"
#include "Dose.h"
#include "DrugType.h"
#include "Global.h"
#include "proteome.h"

using namespace std;


/** A class holding hoshen pkpd drug use info.
 *
 * Each human has an instance for each type of drug present in their blood. */
class HoshenDrug : public Drug {
public:
  /** Create a new instance. */
  HoshenDrug (const DrugType*);
  /** Load an instance from a checkpoint. */
  HoshenDrug (const DrugType*, istream& in);
  void write (ostream& out) const;

  virtual double calculateDrugFactor(const ProteomeInstance* infProteome) const;

protected:
  /** Calculate multiplier to decay a concentration by a duration of time
    *
    * @param time Duration in minutes to decay over */
  virtual double decayFactor (double time);
};

#endif
