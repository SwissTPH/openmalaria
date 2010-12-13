/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#ifndef Hmod_PkPdModel
#define Hmod_PkPdModel

// #include "PkPd/Proteome.h"
#include "util/AgeGroupInterpolation.h"
#include "Global.h"

#include <fstream>
using namespace std;

class UnittestUtil;

namespace OM { namespace PkPd {
    
    using util::AgeGroupInterpolation;

/** Encapsulates both the static operations for PKPD models and the per-human
 * drug proxies.
 * 
 * This is an abstract base class; this helps enforce virtual methods in base
 * classes match (correct name and parameters) since sub-classes not correctly
 * overriding all abstract virtual functions will be abstract and can not be
 * instantiated (i.e. used to create objects).
 * 
 * Calling order within a timestep (see doc for medicate for details):
 *  * getDrugFactor() for each infection
 *  * decayDrugs()
 *  * medicate()
 */
class PkPdModel {
public:
  ///@brief Static functions
  //@{
  static void init ();
  static void cleanup ();
  
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
  
  /** Medicate drugs to an individual, which act on infections the following
   * timesteps, until rendered ineffective by decayDrugs().
   *
   * \param drugAbbrev - The drug abbreviation.
   * \param qty        - the quantity in mg.
   * \param time       - Time in days since start of this time step to medicate at
   * \param ageYears        - Age of human in years
   * 
   * Due to the fact we're using a discrete timestep model, the case-management
   * update (calling medicate) and within-host model update (calling
   * getDrugFactor) cannot [easily] have immediate effects on each other. The
   * implementation we use is that the within-host model update (calculating
   * new infection densities) happens first; hence medicate() will always be
   * called after getDrugFactor in a timestep, and a time of zero means the
   * dose has effect from the start of the following timestep. */
  virtual void medicate(string drugAbbrev, double qty, double time, double ageYears) =0;
  /** Medicate via IV. Mostly as for medicate(). End-time of IV period is passed
   * (time at which concentration is added to use oral effect calculation code).
   */
  virtual void medicateIV(string drugAbbrev, double qty, double duration, double endTime) =0;
  
  /// Called each timestep immediately after the drug acts on any infections.
  virtual void decayDrugs () =0;
  
  /** This is how drugs act on infections.
   *
   * Each timestep, on each infection, the parasite density is multiplied by
   * the return value of this infection. The WithinHostModels are responsible
   * for clearing infections once the parasite density is negligible. */
  virtual double getDrugFactor (uint32_t proteome_ID) =0;
  
  virtual uint32_t new_proteome_ID () =0;
  
  enum ActiveModel {
      NON_PKPD = 0,
      HOSHEN_PKPD,
      LSTM_PKPD
  };
  
protected:
  virtual void checkpoint (istream& stream) =0;
  virtual void checkpoint (ostream& stream) =0;
  
  /** Weight model. Currently looks up a weight dependant on age from a table
   * in an entirely deterministic way.
   *
   * @param ageGroupData Age group for weight data
   * @param ageYears Age in years
   * @param hetMult Multiplies age to introduce heterogeneity
   * @returns Mass in kg */
  static inline double ageToWeight (double ageYears, double hetMult) {
      return (*weight)( ageYears ) * hetMult;
  }
  
  static double hetWeightMultStdDev;
  
protected:
    static double minHetWeightMult;
private:
    static AgeGroupInterpolation* weight;
    
    /// Which model is in use (set by init())
    static ActiveModel activeModel;
    
    friend class ::UnittestUtil;
};

} }
#endif