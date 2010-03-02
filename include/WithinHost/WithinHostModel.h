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

#ifndef Hmod_WithinHost_Model
#define Hmod_WithinHost_Model

#include "Global.h"
#include "WithinHost/Infection.h"
#include "Survey.h"

#include <list>

using namespace std;

namespace OM { namespace WithinHost {
    
/*! Within Host Model abstract class.
 * Dont forget to create friend << and >> for subclasses.
 */
class WithinHostModel {
public:
  /// @brief Static methods
  //@{
  /// Initialise static parameters
  static void init();
  
  /// Free memory
  static void clear();
  
  /// Create an instance using the appropriate model
  static WithinHostModel* createWithinHostModel ();
  //@}
  
  /// @brief Constructors, destructors and checkpointing functions
  //@{
  WithinHostModel();
  virtual ~WithinHostModel() {}
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      checkpoint (stream);
  }
  //@}
  
  void summarize(Survey& survey, SurveyAgeGroup ageGroup);
  
  /// Create a new infection within this human
  virtual void newInfection() =0;
  /** Conditionally clears all infections. Not used with the PK/PD model.
   *
   * If IPT isn't present, it just calls clearAllInfections(); otherwise it
   * uses IPT code to determine whether to clear all infections or do nothing
   * (isSevere is only used in the IPT case). */
  virtual void clearInfections (bool isSevere);
  
  /** Medicate drugs (wraps drug's medicate).
   *
   * @param drugAbbrev	abbrevation of drug name (e.g. CQ, MF)
   * @param qty		Quantity of drug to administer in mg
   * @param time		Time relative to beginning of timestep to medicate at, in days (less than 1 day)
   * @param age		Age of human in years
   */
  virtual void medicate(string drugAbbrev, double qty, double time, double age) {}

  /** Update the parasite densities of infections.
   *
   * @param ageInYears Age of human
   * @param BSVEfficacy Efficacy of blood-stage vaccine */
  virtual void calculateDensities(double ageInYears, double BSVEfficacy) =0;
  
  bool parasiteDensityDetectible() const {
    return totalDensity > detectionLimit;
  }
  
  inline double getTotalDensity() const {return totalDensity;}
  inline double getTimeStepMaxDensity() const {return timeStepMaxDensity;}
  
  ///@brief Only do anything when IPT is present:
  //@{
  /// Conditionally set last SP dose
  virtual void IPTSetLastSPDose (int agetstep, SurveyAgeGroup ageGroup) {}
  /// Prescribe IPTi with probability compliance. Only called if IPT present.
  virtual void IPTiTreatment (SurveyAgeGroup ageGroup);
  //@}
  
  ///@brief Immunity model
  //@{
  /// Called to effect some penalty on immunity − but what? Please document.
  void immunityPenalisation();
  
protected:
  /** Updates for the immunity model − assumes _cumulativeh and _cumulativeY
   * have already been incremented.
   * 
   * Applies decay of immunity against asexual blood stages, if present. */
  void updateImmuneStatus();
  
  /** For summarizing:
   * @returns Total number of infections.
   * @param patentInfections Out param: the number of patent infections
	    (only set if return-value is non-zero). */
  virtual int countInfections (int& patentInfections) =0;
  
  virtual void checkpoint (istream& stream);
  virtual void checkpoint (ostream& stream);
  
  //!innate ability to control parasite densities
  double _innateImmSurvFact;
  
  /** Number of infections received since birth.
   *
   * In the Empirical and Dummy WH models this is updated instantly; in the
   * Descriptive model the update doesn't take effect until 2 timesteps after
   * infection (avoids changing results). TODO: confirm intended behaviour. */
  double _cumulativeh;
  //!Cumulative parasite density since birth
  double _cumulativeY;
  //!cumulativeY from previous timestep
  double _cumulativeYlag;
  //@}
  
  /** Literally just removes all infections in an individual.
   *
   * Normally clearInfections() would be called instead, which, when IPT is not
   * active, just calls this function (although this needs to be changed for
   * PK_PD integration). */
  virtual void clearAllInfections() =0;
  
  //!multiplicity of infection
  int _MOI;
  
  /// Total asexual blood stage density (sum of density of infections).
  double totalDensity;
  
  /** Maximum parasite density of any infection during the previous interval.
   *
   * With 5-day timesteps, this is not just the maximum density of any infection
   * at the end of the timestep, but something designed to emulate the maximum
   * of 5 daily samples. */
  double timeStepMaxDensity;
  
  ///@brief Static parameters, set by init()
  //@{
//Standard dev innate immunity for densities
  static double sigma_i;
// contribution of parasite densities to acquired immunity in the presence of fever
  static double immPenalty_22;
/*
  Remaining immunity against asexual parasites(after time step, each of 2 components y and h)
  This variable decays the effectors cumulativeH and cumulativeY in a way that their
  effects on densities (1-Dh and 1-Dy) decay exponentially.
*/
  static double asexImmRemain;
/*
  Remaining immunity against asexual parasites(after each time step, each of 2 components y and h)
  This variable decays the effectors cumulativeH and cumulativeY exponentially.
*/
  static double immEffectorRemain;
  
  /*
  The detection limit (in parasites/ul) is currently the same for PCR and for microscopy
  TODO: in fact the detection limit in Garki should be the same as the PCR detection limit
  The density bias allows the detection limit for microscopy to be higher for other sites
  */
  static double detectionLimit;
  
  /** The maximum number of infections a human can have. The only real reason
   * for this limit is to prevent incase bad input from causing the number of
   * infections to baloon stupidly. */
  static const int MAX_INFECTIONS = 21;
  //@}
};

} }
#endif
