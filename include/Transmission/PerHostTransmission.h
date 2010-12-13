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
#ifndef Hmod_PerHostTransmission
#define Hmod_PerHostTransmission

#include "util/random.h"
#include "Transmission/Vector/HostCategoryAnopheles.h"
#include "util/AgeGroupInterpolation.h"

namespace OM { namespace Transmission {
    
class HostMosquitoInteraction;
class TransmissionModel;
using util::AgeGroupInterpolation;

/** Contains TransmissionModel parameters which need to be stored per host.
 *
 * Currently many members are public and directly accessed. */
class PerHostTransmission
{
public:
  /// @brief Static member functions
  //@{
  /** Static initialisation. */
  static void init ();
  /** Static cleanup. */
  static void cleanup ();
  
  /** Calculates the adjustment for body size in exposure to mosquitoes,
   * relative to an average adult.
   * 
   * The bites are assumed proportional to average surface area for hosts of
   * the given age. Linear interpolation is used to calculate this from the
   * input array of surface areas. 
   * 
   * @param ageYears Age of host
   * @return the ratio of bites received by the host to the average for an adult 
   *
   * This is the age factor of availiability; mean output should be
   *  1.0/ageCorrectionFactor.
   * 
   * Also has a switch to put individuals entirely outside transmission. */
  inline double relativeAvailabilityAge (double ageYears) const {
    return outsideTransmission ? 0.0 :
	(*relAvailAge)( ageYears );
  }
  //@}
  
  ///@brief Initialisation / checkpionting
  //@{
  PerHostTransmission ();
  void initialise (TransmissionModel& tm, double availabilityFactor);
  //@}
  
  /** @brief Get model parameters for species[speciesIndex].
   *
   * @param speciesIndex = Index in species list of this mosquito type. */
  //@{
  /** Convenience version of entoAvailabilityPartial()*getRelativeAvailability()
   *
   * Mean should be same as entoAvailabilityHetVecItv(). */
    inline double entoAvailabilityFull (
	const HostCategoryAnopheles& base,
	size_t speciesIndex,
	double ageYears,
	double ageCorrectionFactor
    ) const {
	return entoAvailabilityHetVecItv (base, speciesIndex)
	    * relativeAvailabilityAge (ageYears)
	    * ageCorrectionFactor;
  }
  /** @brief Availability of host to mosquitoes (α_i) excluding age factor.
   *
   * (Includes heterogeneity, intervention, and human-to-vector availability
   * rate factors.)
   * 
   * Assume mean is human-to-vector availability rate factor. */
  double entoAvailabilityHetVecItv (const HostCategoryAnopheles& base, size_t speciesIndex) const;
  /** Probability of a mosquito succesfully biting a host (P_B_i). */
  double probMosqBiting (const HostCategoryAnopheles& base, size_t speciesIndex) const;
  /** Probability of a mosquito succesfully finding a resting
   * place after biting and then resting (P_C_i * P_D_i). */
  double probMosqResting (const HostCategoryAnopheles& base, size_t speciesIndex) const;
  //@}
  
  /** Get the availability of this host to mosquitoes relative to an average
   * adult (including heterogeneity and age effects).
   *
   * Used to drive a simulation from an input EIR.
   * Is relativeAvailabilityHet()*relativeAvailabilityAge(ageYears).
   * 
   * Mean output is less than 1.0 (roughly 1.0/ageCorrectionFactor).
   */
  inline double relativeAvailabilityHetAge (double ageYears) const {
    return _relativeAvailabilityHet
        * relativeAvailabilityAge (ageYears);
  }
  /** Relative availability of host to mosquitoes excluding age factor.
   *
   * (ONLY for HeterogeneityWorkaroundII, and documentation purposes.)
   * Assume mean is 1.0. */
  inline double relativeAvailabilityHet () const {
    return _relativeAvailabilityHet;
  }
  
  /** Set true to remove human from transmission. Must set back to false
   * to restore transmission. */
  inline void removeFromTransmission (bool s){
      outsideTransmission = s;
  }
  
  /// Give individual a new ITN as of time timeStep.
  inline void setupITN () {
    timestepITN = Global::simulationTime;
  }
  /// Give individual a new IRS as of time timeStep.
  inline void setupIRS () {
    timestepIRS = Global::simulationTime;
  }
  /// Give individual a new VA intervention as of time timeStep.
  inline void setupVA () {
    timestepVA = Global::simulationTime;
  }
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      species & stream;
      _relativeAvailabilityHet & stream;
      outsideTransmission & stream;
      timestepITN & stream;
      timestepIRS & stream;
      timestepVA & stream;
  }
  
  
private:
  vector<HostMosquitoInteraction> species;
  
  // Heterogeneity factor in availability; this is already multiplied into the
  // entoAvailability param stored in HostMosquitoInteraction.
  double _relativeAvailabilityHet;
  
  // Determines whether human is outside transmission
  bool outsideTransmission;
  
  // (simulationTime - timestepXXX) is the age of the intervention.
  // timestepXXX = TIMESTEP_NEVER means intervention has not been deployed.
  
  int timestepITN;
  int timestepIRS;
  int timestepVA;
  
  static AgeGroupInterpolation* relAvailAge;
};

/** Data needed for each human which is per-mosquito species. */
class HostMosquitoInteraction
{
  friend class PerHostTransmission;
  
public:
  /** In lieu of a constructor initialises elements, using the passed base to
   * get baseline parameters. */
  void initialise (HostCategoryAnopheles& base, double availabilityFactor);
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      entoAvailability & stream;
      probMosqBiting & stream;
      probMosqFindRestSite & stream;
      probMosqSurvivalResting & stream;
  }
  
private:
  ///@brief Rate/probabilities before interventions. See functions.
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
};

} }
#endif
