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

#ifndef Hmod_LSTMDrugType
#define Hmod_LSTMDrugType

#include <string>
#include <deque>
#include <map>
#include <vector>
#include "Dose.h"
#include "Global.h"
#include "PkPd/Proteome.h"
#include "DrugType.h"

using namespace std;

namespace OM { namespace PkPd {
    
    /** Per drug, per age, etc. drug parameters. */						// Add in my parameters and their documentation here
    struct LSTMDrugParameters {
	 /*PD parameters required*/
	double max_killing_rate;			/// Maximal drug killing rate per day
	double IC50;						/// Concentration with 50% of the maximal parasite killing
	double slope;						/// Slope of the dose response curve
	/*PK parameters required*/
	double elimination_rate_constant;	/// Terminal elimination rate constant. Found using ln(2)/half_life
	double vol_dist;					/// Volume of distribution (l/kg)
    };
    
/** Information about each (type of) drug (rather than each use of a drug).
 *
 * Static data contains a list of all available drug types.
 * 
 * No DrugType data is checkpointed, because it is loaded by init() from XML
 * data. (Although if it cannot be reproduced by reloading it should be
 * checkpointed.) */
class LSTMDrugType : public DrugType {
public:
  ///@brief Static functions
  //@{
  /** Initialise the drug model. Called at start of simulation. */
  //TODO: data from XML.
  static void init ();
  
  
  ///@brief Non static (per instance) functions
  //@{
  /** Create a new DrugType.
   *
   * @param name	Name of the drug
   * @param abbreviation	Abbreviated name (e.g. CQ)
   */
  LSTMDrugType (string name, string abbreviation);
  ~LSTMDrugType ();
  
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
    LSTMDrugParameters parameters;
  
  // Allow the Drug class to access private members
  friend class Drug;
  friend class LSTMDrug;
};

} }
#endif