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

#ifndef Hmod_DrugType
#define Hmod_DrugType

#include <string>
#include <deque>
#include <map>
#include <vector>
#include "Dose.h"
#include "Global.h"
#include "PkPd/Proteome.h"

using namespace std;

namespace OM { namespace PkPd {
    
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
  //TODO: data from XML.
  static void init ();
  
  //! Adds a new drug type to the list
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
   */
  DrugType (string name, string abbreviation);
  ~DrugType ();
  //@}
  
protected:
  // The list of available drugs. Not checkpointed; should be set up by init().
  static map<string,DrugType*> available;
  
  //BEGIN Drug-type fields (same for all drugs of same type)
  //! The drug abbreviated name, used for registry lookups.
  string abbreviation;
  //! The drug name.
  string name; 
  
  // Allow the Drug class to access private members
  friend class Drug;
  friend class HoshenDrug;
  friend class IhKwDrug;
};

} }
#endif