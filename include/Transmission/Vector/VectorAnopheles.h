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
#include "Survey.h"
#include "Transmission/Vector/HostCategoryAnopheles.h"
#include "Transmission/Vector/HostCategoryAnophelesNonHumans.h"
#include "Transmission/Vector/HostCategoryAnophelesHumans.h"
#include "Transmission/PerHostTransmission.h"
#include <list>
#include <vector>

namespace OM {
    namespace Host {
	class Human;
    }
namespace Transmission {

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
 * the VectorTransmission class.
 * 
 * Variable names largely come from Nakul Chitnis's paper:
 * "A mathematical model for the dynamics of malaria in mosquitoes feeding on
 * a heterogeneous host population" (3rd Oct. 2007). */

using namespace std;

class VectorAnopheles
{
public:
#ifdef OMV_CSV_REPORTING
    VectorAnopheles (const TransmissionModel* const tm,ofstream& csvR) :
	transmissionModel(tm),
	partialEIR(0.0),
	larvicidingEndStep (std::numeric_limits<int>::max()),
	larvicidingIneffectiveness (1.0),
	timestep_N_v0(0.0), timestep_N_v(0.0), timestep_O_v(0.0), timestep_S_v(0.0),
	csvReporting(&csvR)
    {}
#else
    VectorAnopheles (const TransmissionModel* const tm) :
	transmissionModel(tm),
	partialEIR(0.0),
	larvicidingEndStep (std::numeric_limits<int>::max()),
	larvicidingIneffectiveness (1.0),
	timestep_N_v0(0.0), timestep_N_v(0.0), timestep_O_v(0.0), timestep_S_v(0.0)
    {}
#endif
  
  ///@brief Initialisation and destruction
  //@{
  /** Called to initialise variables instead of a constructor.
   *
   * @param anoph Data structure from XML to use
   * @param sIndex Index in VectorTransmission.species of this class.
   * @param EIR In/out parameter: the EIR used for the pre-intervention phase. Units: innoculations.
   */
  string initialise (const scnXml::Anopheles& anoph, size_t sIndex, vector<double>& EIR, map<string, double>& nonHumanHostsPopulations, int populationSize);
  


  /** Initialise a few more variables (mosqEmergeRate, forcedS_v), which depend
   * on the human population structure (when not loading from a checkpoint).
   * 
   * @param sIndex Index in VectorTransmission.species of this class.
   * @param population The human population
   * @param populationSize Number of humans (use instead of population.size())
   *
   * Can only usefully run its calculations when not checkpointing, due to
   * population not being the same when loaded from a checkpoint. */
  void setupNv0 (size_t sIndex, const std::list<Host::Human>& population, int populationSize);
  
  /** Called to free memory instead of a destructor. */
  void destroy ();
  
  /** Work out whether another interation is needed for initialisation and if
   * so, make necessary changes.
   * 
   * @returns true if another iteration is needed. */
  bool vectorInitIterate ();
  
  /** Calls calMosqEmergeRate() and initialises arrays. */
  //void initMainSimulation (size_t sIndex, const std::list<Host::Human>& population, int populationSize, vector<double>& kappa);
  //@}
  
  /** Called per time-step. Does most of calculation of EIR.
   *
   * @param population The human population; so we can sum up availability and
   *	infectiousness.
   * @param simulationTime
   * @param sIndex Index of the type of mosquito in per-type/species lists.
   * @param isDynamic True to use full model; false to drive model from current contents of S_v. */
  void advancePeriod (const std::list<Host::Human>& population, int simulationTime, size_t sIndex, bool isDynamic);
  
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
  /** Check all required intervention descriptions are present. */
  inline void checkInterventionDescriptions (string species) {
      humanBase.checkInterventionDescriptions (species);
  }
  
  void intervLarviciding (const scnXml::LarvicidingAnopheles&);
  
  /// Write some per-species summary information.
  void summarize (const string speciesName, Survey& survey);
  
  /// Checkpointing
  //Note: below comments about what does and doesn't need checkpointing are ignored here.
  template<class S>
  void operator& (S& stream) {
    humanBase & stream;
    mosqSeekingDeathRate & stream;	
    mosqSeekingDuration & stream;
    mosqRestDuration & stream;
    EIPDuration & stream;
    probMosqSurvivalOvipositing & stream;
    nonHumanHosts & stream;
    EIRRotateAngle & stream;
    FSRotateAngle & stream;
    FSCoeffic & stream;
    mosqEmergeRate & stream;
    forcedS_v & stream;
    annualS_v & stream;
    sumAnnualForcedS_v & stream;
    initNv0FromSv & stream;
    initNvFromSv & stream;
    N_v_length & stream;
    P_A & stream;
    P_df & stream;
    P_dif & stream;
    N_v & stream;
    O_v & stream;
    S_v & stream;
    fArray & stream;
    ftauArray & stream;
    partialEIR & stream;
    larvicidingEndStep & stream;
    larvicidingIneffectiveness & stream;
    timestep_N_v0 & stream;
    timestep_N_v & stream;
    timestep_O_v & stream;
    timestep_S_v & stream;
    initP_A & stream;
    P_A1 & stream;
    P_An & stream;
    mosqLaidEggsSameDayProp & stream;
    probMosqSurvivalFeedingCycle & stream;
  }
  
private:
    /** Reference back to TransmissionModel base. */
    const TransmissionModel* transmissionModel;
    
  /** @brief Baseline parameters which may be varied per host
   *
   * Includes model parameters which may be varied per-individual to account
   * for interventions and innate resistances, and intervention effect
   * descriptions.
   * 
   * Read from XML by initialise; no need to checkpoint. */
  //@{
  HostCategoryAnophelesHumans humanBase;
  //@}
  
  double initP_A;
  double P_A1;
  double P_An;

  /** sets the PAS (See Document "Parameter Values for Transmission model" (Chitnis, Smith and Schapira, 4.3.2010))
   *  initPA : Probability that a mosquito does not find a host and does not die in one night of searching
   *  P_A1 : Probability that a mosquito encounters a human on a given night.
   *  P_An : Probability that a mosquito encounters a non human host on a given night.
   *
   */
  void setPAS();

  /** returns the human ento availability, calculated from PA, PA1, mosqSeekingDuration and population size.
   *
   *  In previous versions the ento availability was to be explicitly given in the scenario. Because of the
   *  dependence of the ento availability value with the population size, we had to change the ento availability
   *  whenever we changed the population size.
   *
   * @param human population size;
   * @return human ento availability;
   *
   */
  double getHumanEntoAvailability(int populationSize);

    /** returns the non human ento availability for a given type of non human host,
     *  calculated from init_PA, P_An, mosqSeekingDuration, population size and relativeEntoAvailability.
     *
     *  In previous versions the ento availability was to be explicitly given in the scenario. Because of the
     *  dependence of the ento availability value with the population size, we had to change the ento availability
     *  whenever we changed the population size.
     *
     *  If only one type of non human host is given in the scenario, then relativeEntoAvailability = 1.
     *
     * @param non human host population size;
     * @return non human host ento availability;
     *
     */
  double getNonHumanEntoAvailability(int populationSize, double relativeEntoAvailability);

  /** Proportion of host-seeking parous mosquitoes that have laid eggs same day*/
  double mosqLaidEggsSameDayProp;

  /** Probability that a mosquito survives a feeding cycle */
  double probMosqSurvivalFeedingCycle;

  /** @brief Parameters which may vary per mosquito species
   *
   * Read from XML by initialise; no need to checkpoint. */
  //@{
  /** Death rate of mosquitoes while host-seeking (μ_vA).
   *
   * TODO: the model could be extended to allow this and mosqSeekingDuration to vary over the year.
   * 
   * Unit: animals/day. */
  double mosqSeekingDeathRate;
  
  /** Duration of host-seeking per day; the maximum fraction of a day that a
   * mosquito would spend seeking (θ_d). */
  double mosqSeekingDuration;
  
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
   * When calcFourierEIR is used to produce an EIR from this over 365
   * (DAYS_IN_YEAR) elements, the resulting EIR has units innoculations.
   * 
   * fcEir must have odd length and is ordered: [a0, a1, b1, ..., an, bn].
   * FSCoeffic[0] needs checkpointing, the rest doesn't. */
  vector<double> FSCoeffic;
  
  /** Emergence rate of new mosquitoes, for every day of the year (N_v0).
   * Units: Animals per day. Length: Global::DAYS_IN_YEAR.
   * 
   * Should be checkpointed. */
  vector<double> mosqEmergeRate;
  
  /* Parameters and partial (derived) parameters from model */
  
  /** @brief S_v used to force an EIR during vector init.
   * Length: Global::DAYS_IN_YEAR
   * 
   * Should be checkpointed. */
  vector<double> forcedS_v;
  
  /** Used by vectorInitIterate to calculate scaling factor.
  *
  * Length of annualS_v is Global::DAYS_IN_YEAR. Checkpoint.
  * 
  * Units of both should be innoculations. */
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
   * Length (fArray): EIPDuration - mosqRestDuration + 1 (θ_s - τ + 1)
   * Length (ftauArray): EIPDuration (θ_s)
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
  
  /** Variables tracking data to be reported. */
  double timestep_N_v0, timestep_N_v, timestep_O_v, timestep_S_v;

  /** map that has populations size as value and non human host type name as key */
  map<string, double> nonHumansHostsPopulations;


  /* Functions */
  
  /** This subroutine converts ShortArray to a vector<double> of length
   * Global::DAYS_IN_YEAR by copying and duplicating elements to fill the gaps. */
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
  //TODO: replace this reporting system
  //NOTE: this is not checkpointed
  ofstream* csvReporting;
#endif
};

} }
#endif
