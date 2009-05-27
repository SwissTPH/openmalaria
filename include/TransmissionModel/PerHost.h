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
  
  PerHostTransmission (TransmissionModel&);
  PerHostTransmission (istream& in, TransmissionModel&);
  void write (ostream& out) const;
  
  void summarize (Summary&, double age);
  
  /*! get the number of infections for a specific human at a given time step 
   *
  1. Calculates h from the EIR measured on adults where 
  h is the expected number of epidemiological inoculations 
  2 Calculates the updated values of the pre-erythrocytic exposure and 
  passes this back to the calling routine 
  Requires the five-day EIR, adjusted for age as input. 
  cumEIR: is the pre-erythrocytic exposure; 
  efficacy: Efficacy of a pre-erythrocytic vaccine 
  \param *cumEIR cumulative EIR (measures pre-erthrocytic immunity) 
  \param efficacy efficacy of a pre-erythrocytic vaccine 
  \param age_adj_EIR Expected number of inoculations adjusted for age of the host 
  \param baseAvailToMos Host-specific availability 
  */ 
  virtual double getExpectedNumberOfInfections (Human& human, double age_adj_EIR);
  
  //! Calculate the number of new infections to introduce via a stochastic process
  int numNewInfections(double expectedInfectionRate, double expectedNumberOfInfections);
  
  //NOTE: may need to be within PerHostTransmission if some intervention parameters are moved here.
  //NOTE: Since only the product of these is usually required, could perhaps be optimised
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
  
  vector<HostMosquitoInteraction> species;
  //@}
  
  /* Shape constant of (Gamma) distribution of availability
  real, parameter :: BaselineAvailabilityGammaShapeParam =1.0 */
  static double BaselineAvailabilityShapeParam;
  
private:
  //NOTE: some intervention data such as timestep of use could be stored here:
  /** simulationTime - dateOfUse is the age of the intervention.
  * 
  * This is the date of last use.
  int timestepITN;
  int timestepIRS; */
  
  
  //VARIABLES INCLUDED IN CORE GETs of number of infections 
  //! The average proportion of bites from sporozoite positive mosquitoes resulting in infection. 
  /*! 
  This is computed as 0.19 (the value S from a neg bin mass action model fitted 
  to Saradidi data, divided by 0.302 (the ratio of body surface area in a 
  0.5-6 year old child (as per Saradidi) to adult) 
  \sa getExpectedNumberOfInfections() 
  */ 
  static const double susceptibility;
  
  //!Steepness of relationship between success of inoculation and Xp in Phase A model 
  /*! 
  \sa getExpectedNumberOfInfections(),Sinf,Simm,Xstar_p,Estar 
  */ 
  static double gamma_p; 
  
  //!Lower limit of success probability of inoculations at high exposure in Phase A model 
  /*! 
  \sa getExpectedNumberOfInfections(),gamma_p,Simm,Xstar_p,Estar 
  */ 
  static double Sinf; 
  
  //!Lower limit of success probability of inoculations in immune individuals in Phase A model 
  /*! 
  \sa getExpectedNumberOfInfections(),gamma_p,Sinf,Xstar_p,Estar 
  */ 
  static double Simm; 
  
  //!Critical value of cumulative number of entomologic inoculations in Phase A model 
  /*! 
  \sa getExpectedNumberOfInfections(),gamma_p,Sinf,Simm,Estar 
  */ 
  static double Xstar_p; 
  
  //!Critical value of EIR in Phase A pre-erythrocytic model 
  /*! 
  \sa getExpectedNumberOfInfections(),gamma_p,Sinf,Simm,Xstar_p 
  */ 
  static double Estar; 
  
  //! Describes the shape of the Infectionrate distribution, related to the baseline availabilty distr. 
  static double InfectionrateShapeParam;
};

/** Data needed for each human which is per-mosquito species. */
class HostMosquitoInteraction
{
  friend class PerHostTransmission;
  
public:
  /** In lieu of a constructor initialises elements, using the passed base to
   * get baseline parameters. */
  void initialise (VectorTransmissionSpecies base);
  
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
