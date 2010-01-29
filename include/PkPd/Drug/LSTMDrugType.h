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
#include "inputData.h"

using namespace std;

namespace OM { namespace PkPd {
    
    /** Per drug, per genotype, PD parameters of drug. */
    struct LSTMDrugPDParameters {
	double max_killing_rate;			/// Maximal drug killing rate per day
	double IC50;						/// Concentration with 50% of the maximal parasite killing
	double slope;						/// Slope of the dose response curve
    };
    
/** Information about each (type of) drug (rather than each use of a drug).
 *
 * Static data contains a list of all available drug types.
 * 
 * No DrugType data is checkpointed, because it is loaded by init() from XML
 * data. (Although if it cannot be reproduced by reloading it should be
 * checkpointed.) */
class LSTMDrugType : public DrugType {	//TODO: check what else is inherited, and possibly remove DrugType base class
public:
  ///@brief Static functions
  //@{
  /** Initialise the drug model. Called at start of simulation. */
  static void init ();
  
  
  ///@brief Non static (per instance) functions
  //@{
  /** Create a new DrugType.
   *
   * @param name	Name of the drug
   * @param abbreviation	Abbreviated name (e.g. CQ)
   */
  LSTMDrugType (const scnXml::Drug& drugData);
  ~LSTMDrugType ();
  //@}
  
private:
    /*PD parameters required - varies with infection genotype*/
    vector<LSTMDrugPDParameters> PD_params;
    
    /*PK parameters required - varies with humans age and severity of disease*/
    double elimination_rate_constant;	/// Terminal elimination rate constant. Found using ln(2)/half_life
    double vol_dist;					/// Volume of distribution (l/kg)
    
  // Allow LSTMDrug to access private members
  friend class LSTMDrug;
  friend inline double drugEffect (const LSTMDrugType& drugType, double& concentration, double duration, double weight_kg, double dose_mg);
};

} }
#endif