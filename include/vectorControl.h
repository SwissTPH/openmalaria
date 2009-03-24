/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef Hmod_VectorControl
#define Hmod_VectorControl
/* This file should contain the headers of all routines that we write in the C
 * program.
 */ 

/* We also include library headers here. */ 
#include "transmissionModel.h"

//#define VectorControl_PRINT_calMosqEmergeRate
//#define VectorControl_PRINT_CalcInitMosqEmergeRate


//! Transmission models, Chitnis et al
class VectorControl : public TransmissionModel {
public:
  VectorControl();
  ~VectorControl();
  
  /** Initialise the main simulation.
   * 
   * Calculates mosquito emergence rate.
   * 
   * \param populationSize The total number of hosts.
   * 
   * Emergence rate calculations assume only one type of host; i.e. it calculates
   * the rate for a stable situation before interventions are introduced. */
  void initMainSimulation (int populationSize); 

  //! get the number of infections for a specific human at a given time step 
  virtual double getExpectedNumberOfInfections (Human& human, double age_adj_EIR);

  /** Calculates EIR (in adults).
   * 
   * \param simulationTime Time since start of simulation . */
  virtual double calculateEIR(int simulationTime, Human& host); 

  /** This needs to be called every interval. */
  virtual void advancePeriod (const std::list<Human>& population, int simulationTime);
  
  //maxDurIntPhaseEIR is the maximum number of days for which an
  //intervention phase EIR is specified
  static const int ifCalcMosqEmergeRate = 0;	// TODO: Move to XML.
  // This should be 1 for the entomological model code to run.
  // But we can set it to 0 to allow the program to run faster by 
  // skipping EntoModel.cpp.

private:
  /* Parameters from model */
  
  /** Emergence rate of new mosquitoes, for every day of the year (N_v0).
   * Units: Animals per day. */
  double mosqEmergeRate[daysInYear];

  /** Death rate of mosquitoes while host-seeking (μ_vA).
   * Unit: animals/day. */
  double mosqSeekingDeathRate;	// FIXME: varies over time
  
  /** Duration of host-seeking per day; the maximum fraction of a day that a
   * mosquito would spend seeking (θ_d). */
  double mosqSeekingDuration;	// FIXME: varies over time
  
  // NOTE: K_vi = P_vi * (host infected) ?
  /** Probability of an infected host infecting a mosquito and that mosquito
   * becoming infective (assuming it survives), per bite.
   * Assumed constant for all hosts. */
  double P_vi;
  
  /** Probability of a mosquito successfully laying eggs after resting (P_E_i).
   * 
   * Currently assumed constant, although NC's non-autonomous model provides
   * an alternative. */
  double probMosqEggLaying;
  
  /** Duration of resting period for mosquito (τ).
   * Units: days. */
  int mosqRestDuration;
  
  /* Partial (derived) parameters from model */
  
  /** Number of days for which data must be stored to calculate N_v, O_v and
   * S_v.
   * 
   * Should equal EIPDuration + mosqRestDuration to allow values up to
   * θ_s + τ - 1 days back, plus current day. */
  int N_v_length;
  
  /** @brief Probability of a mosquito not finding a host one night.
   * NOTE: should length be N_v_length? */
  double *P_A;
  
  /** @brief P_df and P_dif per-day.
   * NOTE: should length be N_v_length?
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
   * Values are recalculated each step (not stored).
   * Length: EIPDuration (θ_s). */
  vector<double> fArray;
  
  /** Per time-step partial calculation of EIR.
   *
   * See comment in advancePeriod() for details of how the EIR is calculated. */
  double partialEIR;
  
  /* Functions */
  
  //! This subroutine converts vectors of length intervalsPerYear to daysInYear. 
  /*! 
    we expect to put the following this into VectorControl class 
 
    Note that ShortArray is assumed to be a pointer to a double array 
    of length intervalsPerYear. We do not explicitly check this. 
 
    For now, we assume that we will use intervalsPerYear and 
    daysInYear as they are defined in global.f. We do not make this 
    subroutine a general subroutine that converts from a given 
    length to another given length. 
 
    Note: We also assume that: 
    daysInYear = interval*intervalsPerYear. 
	 
    \param FullArray an array of doubles of length daysInYear 
    \param ShortArray a pointer to a array of doubles of length intervalsPerYear 
    \sa daysInYear, interval, intervalsPerYear 
  */ 
  void convertLengthToFullYear (double FullArray[daysInYear], double* ShortArray); 
	
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
  void calMosqEmergeRate (int populationSize); 


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
                                double mosqProbOvipositing,
                                double* FHumanInfectivityInitVector,
                                double* FEIRInitVector);
};
#endif
