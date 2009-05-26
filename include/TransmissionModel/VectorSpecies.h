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

#ifndef Hmod_VectorTransmissionSpecies
#define Hmod_VectorTransmissionSpecies

#include "global.h"
#include <list>

namespace scnXml {
  class Anopheles;
}
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
class VectorTransmissionSpecies
{
public:
  ///@brief Initialisation and destruction
  //@{
  /** Called to initialise variables instead of a constructor.
   *
   * @param anoph Data structure from XML to use
   * @param EIR In/out parameter: the EIR used for the pre-intervention phase.
   */
  void initialise (const scnXml::Anopheles& anoph, vector<double>& EIR);
  
  /** Called to free memory instead of a destructor. */
  void destroy ();
  //@}
  
  /*! get mosquito emergence rates 
   *
   * This routine passes the basic entomological parameters (that are already
   * been read, the EIR, and the human infectivity to mosquitoes (all for one
   * type of host) and calculate the mosquito emergence.
   * 
   * \param populationSize
   * 	Number of hosts of each type.
   * 	Units: Animals. 
   * 	$N_i$ in model. Matrix of size $n \times \theta_p$.
   * 	We assume that the size of the one group in initialization is
   * 	fixed over the cycle.
   * 	Mathematically, we require this parameter to be a positive
   * 	real number, so although this will typically be a natural 
   * 	number, it is not restricted to being one. */
  void calMosqEmergeRate (int populationSize, vector<double>& initialKappa); 
  
  /** Called per time-step. Does most of calculation of EIR.
   *
   * @param sIndex Index of the type of mosquito in per-type/species lists. */
  void advancePeriod (const std::list<Human>& population, int simulationTime, size_t sIndex);
  
  ///@brief Parameters which may vary per mosquito species
  //@{
  /** Emergence rate of new mosquitoes, for every day of the year (N_v0).
   * Units: Animals per day. */
  double mosqEmergeRate[daysInYear];
  
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
  //@}
  
  /** @brief Baseline parameters which may be varied per host
   *
   * These may be varied per-human to account for interventions and innate
   * resistances. */
  //@{
  /** Availability rate (α_i) */
  double entoAvailability;
  
  /** Probability of mosquito successfully biting host (P_B_i) */
  double probMosqBiting;
  
  /** Probability of mosquito escaping human and finding a resting site without
   * dying, after biting the human (P_C_i). */
  double probMosqFindRestSite;
  
  /** Probability of mosquito successfully resting after finding a resting site
   * (P_D_i). */
  double probMosqSurvivalResting;
  //@}
  
  /** Per time-step partial calculation of EIR.
  *
  * See comment in advancePeriod() for details of how the EIR is calculated. */
  double partialEIR;
  
private:
  /* Parameters from model */
  /* Partial (derived) parameters from model */
  
  /** Number of days for which data must be stored to calculate N_v, O_v and
   * S_v.
   * 
   * Should equal EIPDuration + mosqRestDuration to allow values up to
   * θ_s + τ - 1 days back, plus current day. */
  int N_v_length;
  
  /** @brief Probability of a mosquito not finding a host one night.
   * NOTE: needs to be N_v_length long? */
  double *P_A;
  
  /** @brief P_df and P_dif per-day.
   * NOTE: needs to be N_v_length long?
   * 
   * P_df is the probability of a mosquito finding a host and completing a
   * feeding cycle without being killed.
   * 
   * P_dif is the probability of a mosquito finding a host, getting infected,
   * and successfully completing a feeding cycle. */
  double *P_df, *P_dif;
  
  /** Number of host-seeking mosquitos each day; respectively: total number,
   * infected (dela, and infective.
   * Index for each day is day % N_v_length.
   * Length: N_v_length (longer than needed for S_v, but simplifyies code) */
  double *N_v, *O_v, *S_v;
  
  /** Used to calculate recursive functions f and f_τ in NDEMD eq 1.6, 1.7.
   * Values are recalculated each step, only first few elements are stored
   * across steps.
   * Length: EIPDuration (θ_s). */
  //@{
  vector<double> fArray;
  vector<double> ftauArray;
  //@}
  
  
  /** @brief Parameters used during the initialisation phase. */
  //@{
  // NOTE: only needs to be stored so calMosqEmergeRate doesn't have to recalculate.
  /** EIR for this species to be used during the initialisation phase.
   *
   * Length: Global::intervalsPerYear */
  vector<double> speciesEIR;
  
  /** The filename to which emergence rates are loaded & saved. */
  string emergenceRateFilename;
  //@}
  
  /* Functions */
  
  /** This subroutine converts ShortArray of length intervalsPerYear to
   * FullArray by copying and duplicating elements to fill the gaps.
   *
   * It also rotates FullArray slightly, since the indecies of some input
   * arrays used to be one-based but are now zero based while the use of the
   * output is unchanged. NOTE: may not be correct, but this preserves the
   * previous situation. */
  void convertLengthToFullYear (double FullArray[daysInYear], vector<double>& ShortArray); 


  /** calcInitMosqEmergeRate() calculates the mosquito emergence rate given
   * all other parameters.
   *
   * We use a periodic version of the model described in "A Mathematical Model 
   * for the Dynamics of Malaria in Mosquitoes Feeding on a Heteregeneous Host
   * Population". The periodic model still needs to be written as a paper. We will
   * change these comments to refer to the approprirate paper when it is ready.
   *
   * The entomological model has a number of input parameters, including the
   * mosquito emergence rate, $N_{v0}$, and a number of output parameters, 
   * including the entomological inoculation rate, $\Xi_i$. The model produces
   * equations for $\Xi_i$ as a function of $N_{v0}$ and the other parameters.
   * However, in this function, we assume that all parameters, except $N_{v0}$ 
   * are known, and $\Xi_i$ is known. We then use these parameters, with $\Xi_i$ 
   * to calculate $N_{v0}$. The equations for $\Xi_i$ are linear in terms of 
   * $N_{v0}$ so there is a unique solution for $N_{v0}$. 
   *
   * This routine first shows the existence of a unique globally asymptotically 
   * stable periodic orbit for the system of equations describing the periodically
   * forced entomological model (for a given set of parameter values, including the
   * mosquito emergence rate). It then compares the number of infectious host-seeking
   * mosquitoes for this periodic orbit to the the number of infectious host-seeking
   * mosquitoes that would result in the given EIR. The routine then iteratively finds
   * the emergence rate that matches the given EIR.
   * 
   * However, we cannot write these equations in the form Ax=b, so we use
   * a root-finding algorithm to calculate $N_{v0}$.
   *
   * This function has a dummy return of 0.
   * 
   * All parameters are IN parameters. */
  double CalcInitMosqEmergeRate(int populationSize,
                                int EIPDuration,
                                int nHostTypesInit, int nMalHostTypesInit,
                                double hostAvailabilityRateInit,
                                double mosqProbBiting,
                                double mosqProbFindRestSite,
                                double mosqProbResting,
                                double* FHumanInfectivityInitVector,
                                double* FEIRInitVector);
  
  /** Given a positive array, originalArray, of length OALength,
   * this routine exponentiates the inverse discrete Fourier 
   * tranform of the first three modes of the natural logarithm of 
   * the array to smooth out the array to produce smoothArray of 
   * length SALength.
   *
   * All elements of originalArray are assumed to be strictly
   * positive.
   *
   * smoothArray is an OUT parameter.
   * originalArray, SALength and OALength are IN parameters.
   * No reason smoothArray and originalArray can't be the same array. */
  void logDFTThreeModeSmooth (double* smoothArray, double* originalArray, int SALength, int OALength); 

  /**
   *  Given a sequence of Fourier coefficients, FC, of odd length,
   *  this routine calculates the exponent of the inverse discrete
   *  Fourier transform into an array, tArray.
   *
   * tArray is an OUT parameter, FC is an IN parameter. */
  void calcInverseDFTExp(vector<double>& tArray, vector<double>& FC);

  /// Shifts elements of rArray clockwise by rAngle.
  void rotateArray(vector<double>& rArray, double rAngle);
};

#endif
