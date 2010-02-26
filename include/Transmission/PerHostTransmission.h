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

#include "WithinHost/WithinHostModel.h"	// for getAgeGroup()
#include "util/random.h"
#include "inputData.h"
#include "Transmission/Vector/HostCategoryAnopheles.h"

namespace OM { namespace Transmission {
    
class HostMosquitoInteraction;
class TransmissionModel;

/** Contains TransmissionModel parameters which need to be stored per host.
 *
 * Currently many members are public and directly accessed. */
class PerHostTransmission
{
public:
  /// @brief Static member functions
  //@{
  /** Static initialisation
   *
   * Pass getInterventions() from inputData. */
  static void initParameters (const scnXml::Interventions&);
  
  //! Calculates the adjustment for body size in exposure to mosquitoes 
  /*! 
  The bites are assumed proportional to average surface area for hosts of the given age. 
  Linear interpolation is used to calculate this from the input array of surface areas. 
  \param ageyrs age in years 
  \return the ratio of bites received by the host to the average for an adult 
   *
   * Age factor of availiability; to be multiplied by partial availability.
   * 
   * Mean output should be 1.0/ageCorrectionFactor. */
  static double relativeAvailabilityAge (double ageyrs) {
    return ageSpecificRelativeAvailability[WithinHost::WithinHostModel::getAgeGroup(ageyrs)];
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
  inline double entoAvailabilityFull (const HostCategoryAnopheles& base, size_t speciesIndex, double ageYears, double ageCorrectionFactor) const {
    return entoAvailabilityHetVecItv (base, speciesIndex)
         * relativeAvailabilityAge (ageYears) * ageCorrectionFactor;
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
  
  /** Get the availability of this host to mosquitoes relative to other hosts.
   * (Includes heterogeneity and age effects; for non-vector model this is the
   * whole availability factor.)
   *
   * Used to drive a simulation from an input EIR.
   * Is relativeAvailabilityHet()*relativeAvailabilityAge(ageYears).
   * 
   * Mean output is less than 1.0 (roughly 1.0/ageCorrectionFactor);
   * this is to avoid changing all NonVector model results. */
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
  
  /** Distribute ITNs to individuals of the correct age (to model ITN
   * distribution along with measles vaccines, etc.). */
  void continousItnDistribution (int ageTSteps);
  
  /// Checkpointing
  template<class S>
  void operator& (S& stream) {
      species & stream;
      _relativeAvailabilityHet & stream;
      timestepITN & stream;
      timestepIRS & stream;
      timestepVA & stream;
      nextItnDistribution & stream;
  }
  
  
private:
  vector<HostMosquitoInteraction> species;
  
  // Heterogeneity factor in availability; this is already multiplied into the
  // entoAvailability param stored in HostMosquitoInteraction.
  double _relativeAvailabilityHet;
  
  // (simulationTime - timestepXXX) is the age of the intervention.
  // timestepXXX = TIMESTEP_NEVER means intervention has not been deployed.
  
  int timestepITN;
  int timestepIRS;
  int timestepVA;
  
  size_t nextItnDistribution;
  
  
  ///@brief Age-group variables for wtprop and ageSpecificRelativeAvailability
  // set by initParameters
  //@{
  /** Average number of bites for each age as a proportion of the maximum.
   *
   * Set by constructor from constants (bsa_prop). */
  static double ageSpecificRelativeAvailability[WithinHost::WithinHostModel::nages];

  //! Proportionate body surface area
 /* 
  The body surface area is expressed as proportions of 0.5*those in 
  the reference age group.In some models we have used calculations of weight and in others surface area, based on 
  Mosteller RD: Simplified Calculation of Body Surface Area. N Engl J Med 1987 Oct 22;317(17):1098 (letter) 
  These values are retained here should they be required for future comparisons 
 */ 
  static const double bsa_prop[WithinHost::WithinHostModel::nages];
  //@}
  
  /// Target ages at which individuals may receive ITNs
  static vector<double> cntItnTargetAgeTStep;
  /// Coverage levels associated with cntItnTargetAgeTStep
  static vector<double> cntItnCoverage;
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
