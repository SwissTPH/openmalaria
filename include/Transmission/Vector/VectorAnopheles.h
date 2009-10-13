/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef Hmod_VectorAnopheles
#define Hmod_VectorAnopheles

#include "Global.h"
#include "Transmission/Vector/HostCategoryAnopheles.h"
#include "Transmission/PerHostTransmission.h"
#include <list>

class Human;

/** Per-species data for vector control.
 *
 * Data in this class is specific to the "species" of anopheles mosquito, where
 * species is used in a relaxed way to mean any variation of anopheles
 * mosquito, not just those types formally recognised as distinct species.
 *
 * A list of this class type is used by the VectorTransmission class to hold
 * (potentially) species-specific per-population data.
 * 
 * Instead of storing static variables in this class, store them in
 * the VectorTransmission class. */
class VectorAnopheles
{
public:
  VectorAnopheles (
#ifdef OMV_CSV_REPORTING
  ofstream& csvR
#endif
  ) :
    partialEIR(0.0),
    larvicidingEndStep (std::numeric_limits<int>::max()),
    larvicidingIneffectiveness (1.0)
#ifdef OMV_CSV_REPORTING
  , csvReporting(&csvR)
#endif
  {}
  
  ///@brief Initialisation and destruction
  //@{
  /** Called to initialise variables instead of a constructor.
   *
   * @param anoph Data structure from XML to use
   * @param sIndex Index in VectorTransmission.species of this class.
   * @param EIR In/out parameter: the EIR used for the pre-intervention phase.
   */
  string initialise (const scnXml::Anopheles& anoph, size_t sIndex, vector<double>& EIR);
  
  void write(ostream& out) const;
  void read(istream& in);
  
  /** Initialise a few more variables (mosqEmergeRate, forcedS_v), which depend
   * on the human population structure.
   * 
   * @param sIndex Index in VectorTransmission.species of this class.
   * @param population The human population
   * @param populationSize Number of humans (use instead of population.size())
   */
  void setupNv0 (size_t sIndex, const std::list<Human>& population, int populationSize);
  
  /** Called to free memory instead of a destructor. */
  void destroy ();
  
  /** Work out whether another interation is needed for initialisation and if
   * so, make necessary changes.
   * 
   * @returns true if another iteration is needed. */
  bool vectorInitIterate ();
  
  /** Calls calMosqEmergeRate() and initialises arrays. */
  //void initMainSimulation (size_t sIndex, const std::list<Human>& population, int populationSize, vector<double>& kappa);
  //@}
  
  /** Called per time-step. Does most of calculation of EIR.
   *
   * @param population The human population; so we can sum up availability and
   *	infectiousness.
   * @param simulationTime
   * @param sIndex Index of the type of mosquito in per-type/species lists.
   * @param isDynamic True to use full model; false to drive model from current contents of S_v. */
  void advancePeriod (const std::list<Human>& population, int simulationTime, size_t sIndex, bool isDynamic);
  
  /** Returns the EIR calculated by advancePeriod().
   * 
   * Could be extended to allow input EIR driven initialisation on a per-species
   * level instead of the whole simulation, but that doesn't appear worth doing.
   * 
   * @param sIndex Index of this in VectorTransmission.species
   * @param host PerHostTransmission of the human requesting this EIR. */
  double calculateEIR (size_t sIndex, PerHostTransmission& host) {
    /* Calculates EIR per individual (hence N_i == 1).
     *
     * See comment in VectorAnopheles.advancePeriod for method. */
    return partialEIR
      * host.entoAvailabilityHetVecItv (humanBase, sIndex)
      * host.probMosqBiting(humanBase, sIndex);	// probability of biting, once commited
  }
  
  inline HostCategoryAnopheles& getHumanBase () {
    return humanBase;
  }
  
  /** Set up intervention descriptions for humans, for this anopheles species.
   *
   * Currently no interventions for non-human hosts, although planned. */
  inline void setInterventionDescription (const scnXml::Anopheles1& intervDesc) {
    humanBase.setInterventionDescription (intervDesc);
  }
  
  void intervLarviciding (const scnXml::LarvicidingAnopheles&);
  
private:
  /** @brief Baseline parameters which may be varied per host
   *
   * Includes model parameters which may be varied per-individual to account
   * for interventions and innate resistances, and intervention effect
   * descriptions.
   * 
   * Read from XML by initialise; no need to checkpoint. */
  //@{
  HostCategoryAnopheles humanBase;
  //@}
  
  /** @brief Parameters which may vary per mosquito species
   *
   * Read from XML by initialise; no need to checkpoint. */
  //@{
  /** Death rate of mosquitoes while host-seeking (μ_vA).
   * Unit: animals/day. */
  double mosqSeekingDeathRate;	// TODO: varies over time
  
  /** Duration of host-seeking per day; the maximum fraction of a day that a
   * mosquito would spend seeking (θ_d). */
  double mosqSeekingDuration;	// TODO: varies over time
  
  /** Duration of resting period for mosquito (τ).
   * Units: days. */
  int mosqRestDuration;
  
  /** Duration of the extrinsic incubation period (sporozoite development time)
  * (θ_s).
  * Units: Days.
  * 
  * Doesn't need checkpointing. */
  int EIPDuration;
  
  /** Probability of a mosquito successfully laying eggs given that it has
   * rested (P_E).
   * 
   * Currently assumed constant, although NC's non-autonomous model provides
   * an alternative. */
  double probMosqSurvivalOvipositing;
  
  /** Non-human host data. Doesn't need checkpointing. */
  NonHumanHostsType nonHumanHosts;
  
  /// Angle (in radians) to rotate series generated by FSCoeffic by, for EIR.
  double EIRRotateAngle;
  //@}
  
  /// Rotation angle (in radians) for emergence rate. Both offset for EIR given in XML file and
  /// offset needed to fit target EIR (delayed from emergence rate). Checkpoint.
  double FSRotateAngle;
  
  /** Fourier coefficients for EIR / forcedS_v series, input from XML file.
   * 
   * Initially used to calculate initialisation EIR, then scaled to calc. S_v.
   * 
   * fcEir must have odd length and is ordered: [a0, a1, b1, ..., an, bn].
   * Doesn't currently need checkpointing. */
  vector<double> FSCoeffic;
  
  /** Emergence rate of new mosquitoes, for every day of the year (N_v0).
   * Units: Animals per day. Length: daysInYear.
   * 
   * Set by setupNv0, then adjusted; should be checkpointed. */
  vector<double> mosqEmergeRate;
  
  /* Parameters and partial (derived) parameters from model */
  
  /** @brief S_v used to force an EIR during vector init.
   * Length: daysInYear
   * 
   * Should be checkpointed. */
  vector<double> forcedS_v;
  
  /** Used by vectorInitIterate to calculate scaling factor.
  *
  * Length of annualS_v is daysInYear. Checkpoint. */
  vector<double> annualS_v;
  double sumAnnualForcedS_v;	///< ditto
  
  /** Conversion factor from forcedS_v to mosqEmergeRate.
   *
   * Also has another temporary use between initialise and setupNv0 calls.
   * Should be checkpointed. */
  double initNv0FromSv;	///< ditto
  
  /** Conversion factor from forcedS_v to (initial values of) N_v.
   * Should be checkpointed. */
  double initNvFromSv;
  
  /** N_v_length-1 is the number of previous days for which some parameters are
   * stored: P_A, P_df, P_dif, N_v, O_v and S_v. This is longer than some of
   * the arrays need to be, but simplifies code with no real impact.
   * 
   * Should equal EIPDuration + mosqRestDuration to allow values up to
   * θ_s + τ - 1 days back, plus current day.
   * 
   * Set by initialise; no need to checkpoint. */
  int N_v_length;
  
  /** @brief Parameter arrays N_v_length long.
   *
   * P_A, P_df, P_dif, N_v, O_v and S_v are set in advancePeriod().
   * 
   * They should be checkpointed. */
  //@{
  /** @brief Probability of a mosquito not finding a host one night. */
  vector<double> P_A;
  
  /** @brief P_df and P_dif per-day.
   * 
   * P_df is the probability of a mosquito finding a host and completing a
   * feeding cycle without being killed.
   * 
   * P_dif is the probability of a mosquito finding a host, getting infected,
   * and successfully completing a feeding cycle.
   * 
   * HOWEVER, if the initialisation phase is driven by an input EIR and not by
   * vector calculations, then during the initialisation phase, P_dif contains
   * the daily kappa values read from XML for validation purposes. */
  vector<double> P_df, P_dif;
  
  /** Number of host-seeking mosquitos each day; respectively: total number,
   * infected, and infective.
   * Index for each day is day % N_v_length. */
  vector<double> N_v, O_v, S_v;
  //@}
  
  /** @brief Used for calculations within advancePeriod. Only saved for optimisation.
   *
   * Used to calculate recursive functions f and f_τ in NDEMD eq 1.6, 1.7.
   * Values are recalculated each step; only fArray[0] and
   * ftauArray[0..mosqRestDuration] are stored across steps for optimisation.
   * 
   * Length: EIPDuration (θ_s).
   * 
   * Don't need to be checkpointed, but some values need to be initialised. */
  //@{
  vector<double> fArray;
  vector<double> ftauArray;
  //@}
  
  
  /** Per time-step partial calculation of EIR.
  *
  * See comment in advancePeriod() for details of how the EIR is calculated.
  * 
  * Doesn't need to be checkpointed (is recalculated each step). */
  double partialEIR;
  
  /** @brief Simple larviciding intervention.
   *
   * Would need to be checkpointed for main simulation; not used during
   * initialisation period (so can be reinitialised). */
  //@{
  /** Timestep at which larviciding effects dissappear. */
  int larvicidingEndStep;
  /** One-minus larviciding effectiveness. I.e. emergence rate is multiplied by
   * this parameter. */
  double larvicidingIneffectiveness;
  //@}
  
  /* Functions */
  
  /** This subroutine converts ShortArray to a vector<double> of length
   * daysInYear by copying and duplicating elements to fill the gaps. */
  static vector<double> convertLengthToFullYear (vector<double>& ShortArray); 
  
  /** Given an input sequence of Fourier coefficients, with odd length,
   * calculate the exponential of the corresponding fourier series.
   * 
   * Note: output is per-interval in tArray. When length is intervalsPerYear,
   * you may want to scale the output by days-per-interval.
   * 
   * @param tArray Array to fill with EIR values. Length should already be set.
   * @param FC Fourier coefficients (a0, a1,b1, a2,b2, ...).
   * @param rAngle Angle to rotate EIR, in radians: [0,2π] */
  static void calcFourierEIR (vector<double>& tArray, vector<double>& FC, double rAngle);

  /// Shifts elements of rArray clockwise by rAngle.
  static void rotateArray(vector<double>& rArray, double rAngle);
  
  friend class VectorEmergenceSuite;
  friend class VectorAnophelesSuite;
  
#ifdef OMV_CSV_REPORTING
  /// This is used to output infectiousness, etc. as a csv file, when included
  ofstream* csvReporting;
#endif
};

#endif
