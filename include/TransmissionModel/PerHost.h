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

#include "EntoIntervention.h"
#include "TransmissionModel/VectorSpecies.h"

class Summary;
class HostMosquitoInteraction;
class TransmissionModel;

/** Contains TransmissionModel parameters which need to be stored per host.
 *
 * Currently many members are public and directly accessed. */
// TODO: optimise for memory
class PerHostTransmission
{
public:
  /// @brief Static member functions
  //@{
  /// Static initialisation
  static void initParameters ();
  //@}
  
  PerHostTransmission ();
  void initialise (TransmissionModel& tm, double availabilityFactor);
  PerHostTransmission (istream& in, TransmissionModel&);
  void write (ostream& out) const;
  
  //NOTE: may need to be within PerHostTransmission if some intervention parameters are moved here.
  //NOTE: Since only the product of these is usually required, could perhaps be optimised
  /** @brief Get model parameters for species[speciesIndex].
   *
   * @param speciesIndex = Index in species list of this mosquito type. */
  //@{
  /** Availability of host to mosquitoes (α_i).
   *
   * The full availability is
   * entoAvailability() * getRelativeAvailability(human->getAgeInYears()). */
  double entoAvailability (size_t speciesIndex) const;
  /** Probability of a mosquito succesfully biting a host (P_B_i). */
  double probMosqBiting (size_t speciesIndex) const;
  /** Probability of a mosquito succesfully finding a resting
  * place after biting (P_C_i). */
  double probMosqFindRestSite (size_t speciesIndex) const;
  /** Probability of a mosquito succesfully resting (P_D_i). */
  double probMosqSurvivalResting (size_t speciesIndex) const;
  //@}
  
  /** Get the availability of this host to mosquitoes (excl. age factor).
   *
   * For NonVector and initialisation phase with Vector. */
  double entoAvailability () const {
    return _entoAvailability;
  }
  
private:
  vector<HostMosquitoInteraction> species;
  
  // Only used in the non-vector model and initialisation phase of the vector model.
  double _entoAvailability;
  
  //NOTE: some intervention data such as timestep of use could be stored here:
  /** simulationTime - dateOfUse is the age of the intervention.
  * 
  * This is the date of last use.
  int timestepITN;
  int timestepIRS; */
};

/** Data needed for each human which is per-mosquito species. */
class HostMosquitoInteraction
{
  friend class PerHostTransmission;
  
public:
  /** In lieu of a constructor initialises elements, using the passed base to
   * get baseline parameters. */
  void initialise (VectorTransmissionSpecies base, double availabilityFactor);
  
  void read (istream& in);
  void write (ostream& out) const;
  
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
  
  /// Intervention: an ITN (active if netEffectiveness > 0)
  EntoInterventionITN entoInterventionITN;
  /// Intervention: IRS (active if insecticide != 0)
  EntoInterventionIRS entoInterventionIRS;
};

#endif
