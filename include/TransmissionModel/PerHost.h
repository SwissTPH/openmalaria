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

class Summary;

/** Contains TransmissionModel parameters which need to be stored per host.
 *
 * Currently many members are public and directly accessed.
 * 
 * 
 */
class PerHostTransmission
{
public:
  /// @brief Static member functions
  //@{
  /// Static initialisation
  static void initParameters ();
  //@}
  
  PerHostTransmission ();
  
  
  void read (istream& in);
  void write (ostream& out) const;
  
  void summarize (Summary&, double age);
  
  //! Calculate the number of new infections to introduce via a stochastic process
  int numNewInfections(double expectedInfectionRate, double expectedNumberOfInfections);
  
  ///@brief Get model parameters for species[speciesIndex].
  //@{
  /** Availability of host to mosquitoes (α_i). */
  double entoAvailability (size_t speciesIndex) const;
  /** Probability of a mosquito succesfully biting a host (P_B_i). */
  double probMosqBiting (size_t speciesIndex) const;
  /** Probability of a mosquito succesfully finding a resting
  * place after biting (P_C_i). */
  double probMosqFindRestSite (size_t speciesIndex) const;
  /** Probability of a mosquito succesfully resting (P_D_i). */
  double probMosqSurvivalResting (size_t speciesIndex) const;
  //@}
  
  /// @brief Public data members
  //@{
  //!Number of infective bites since birth
  double _cumulativeEIRa;
  //!pinfected: probability of infection (cumulative or reset to zero in massTreatment). Appears to be used only for calculating expected inoculations for the analysis
  //!of pre-erythrocytic immunity.
  double _pinfected;
  //!Baseline availability to mosquitoes
  double _BaselineAvailabilityToMosquitoes;
  
  vector<PerHostPerSpecies> species;
  //@}
  
  /* Shape constant of (Gamma) distribution of availability
  real, parameter :: BaselineAvailabilityGammaShapeParam =1.0 */
  static double BaselineAvailabilityShapeParam;
  
private:
  //NOTE: include intervention data here?
  /** simulationTime - dateOfUse is the age of the intervention.
  * 
  * This is the date of last use. */
  int timestepITN;
  int timestepIRS;
  
  /** Insecticide used. */
  int insecticideITN;
  int insecticideIRS;
  
  //NOTE: or in sub-containers
  /// Intervention: an ITN (active if netEffectiveness > 0)
  EntoInterventionITN entoInterventionITN;
  /// Intervention: IRS (active if insecticide != 0)
  EntoInterventionIRS entoInterventionIRS;
};

/** Data needed for each human which is per-mosquito species. */
class PerHostPerSpecies
{
public:
  PerHostPerSpecies();
  
  friend class PerHostTransmission;
  
private:
  /// Rate/probabilities before interventions. See functions.
  double entoAvailability;
  double probMosqBiting;
  double probMosqFindRestSite;
  double probMosqSurvivalResting;
}

#endif
