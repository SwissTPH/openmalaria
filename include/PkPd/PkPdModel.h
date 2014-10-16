/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef Hmod_PkPdModel
#define Hmod_PkPdModel

// #include "PkPd/Proteome.h"
#include "Global.h"

#include <fstream>
using namespace std;

class UnittestUtil;

namespace scnXml{ class Scenario; }
namespace OM { namespace PkPd {
    
/** Encapsulates both the static operations for PKPD models and the per-human
 * drug proxies.
 * 
 * This is an abstract base class; this helps enforce virtual methods in base
 * classes match (correct name and parameters) since sub-classes not correctly
 * overriding all abstract virtual functions will be abstract and can not be
 * instantiated (i.e. used to create objects).
 * 
 * Calling order each day:
 *  * prescribe()
 *  * medicate()
 *  * getDrugFactor() for each infection
 *  * decayDrugs()
 */
class PkPdModel {
public:
  ///@brief Static functions
  //@{
  static void init ( const scnXml::Scenario& scenario );
  
  // checkpointing of static data: not required since all data is set up by init
  static void staticCheckpoint (istream& stream) {}
  static void staticCheckpoint (ostream& stream) {}
  
  /** Factory function to create a drug interface, type dependant on run-time
   * options.
   * 
   * Currently may return one of: PkPdModel, IhKwPkPdModel. */
  static PkPdModel* createPkPdModel ();
  //@}
  
  /// @brief Constructors, destructor and checkpointing
  //@{
  /// Create a new instance
  PkPdModel () {}
  /// Destroy an instance
  virtual ~PkPdModel () {}
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      checkpoint (stream);
  }
  //@}
  
    /** Feature added for a drug monitoring HACK. Could be replaced with better
     * monitoring code.
     * 
     * Gives the drug concentrations at the start of the day (before medication,
     * where drugs are medicated at hour 0).
     * 
     * @concentrations Table; entries of the form ("LU",0.1) are set for all
     *  modeled drugs (i.e. none unless drugs were recently medicated).
     */
    virtual void getConcentrations(map<string,double>& concentrations) const =0;
  
    /** Prescribe a patient a course of drugs, via the Pk/Pd model
     * 
     * Note: doses sizes are modified according to age via the dosage
     * table given at the time this function is called.
     *
     * @param schedule Index of a treatment schedule
     * @param dosages Index of a dosage table
     * @param age Age of human at start of time step in years
     */
    virtual void prescribe(size_t schedule, size_t dosages, double age) =0;
    
    /** Medicate drugs: human takes prescribed drugs which are to be taken this
     * day.
     * 
     * @param age Age of human in years
     * 
     * Note: poor adherence on the part of the patient is not modeled here; to
     * model, prescribe with a "poor adherence" schedule.
     */
    virtual void medicate(double age) =0;
    
    /** This is how drugs act on infections.
     *
     * Each time step, on each infection, the parasite density is multiplied by
     * the return value of this infection. The WithinHostModels are responsible
     * for clearing infections once the parasite density is negligible. */
    virtual double getDrugFactor (uint32_t proteome_ID) =0;
    
    /** After any resident infections have been reduced by getDrugFactor(),
     * this function is called to update drug levels to their effective level
     * at the end of the day, as well as clear data once drug concentrations
     * become negligible. */
    virtual void decayDrugs () =0;
  
  virtual uint32_t new_proteome_ID () =0;
  
  enum ActiveModel {
      NON_PKPD = 0,
//       HOSHEN_PKPD,   note: this code is no longer maintained or enabled
      LSTM_PKPD
  };
  
protected:
  virtual void checkpoint (istream& stream) =0;
  virtual void checkpoint (ostream& stream) =0;
  
private:
    /// Which model is in use (set by init())
    static ActiveModel activeModel;
    
    friend class ::UnittestUtil;
};

} }
#endif