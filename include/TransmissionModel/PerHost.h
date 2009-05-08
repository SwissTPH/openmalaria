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

/** For now, this is just a data container.
 *
 * Currently members are public; later make private and use friend classes?
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
  
  /** Availability of host to mosquitoes (α_i). */
  double entoAvailability () const;
  /** Probability of a mosquito succesfully biting a host (P_B_i). */
  double probMosqSurvivalBiting () const;
  /** Probability of a mosquito succesfully finding a resting
  * place and resting (P_C_i * P_D_i). */
  double probMosqSurvivalResting () const;
  
  /// @brief Public data members
  //@{
  //!Number of infective bites since birth
  double _cumulativeEIRa;
  //!pinfected: probability of infection (cumulative or reset to zero in massTreatment). Appears to be used only for calculating expected inoculations for the analysis
  //!of pre-erythrocytic immunity.
  double _pinfected;
  //!Baseline availability to mosquitoes
  double _BaselineAvailabilityToMosquitoes;
  
  /// Rate/probabilities before interventions. See functions.
  double _entoAvailability;
  double _probMosqSurvivalBiting;
  double _probMosqSurvivalResting;
  
  /// Intervention: an ITN (active if netEffectiveness > 0)
  EntoInterventionITN entoInterventionITN;
  /// Intervention: IRS (active if insecticide != 0)
  EntoInterventionIRS entoInterventionIRS;
  //@}
  
  /* Shape constant of (Gamma) distribution of availability
  real, parameter :: BaselineAvailabilityGammaShapeParam =1.0 */
  static double BaselineAvailabilityShapeParam;
  
private:
  /* Static private */
  
  static double baseEntoAvailability;
  static double baseProbMosqSurvivalBiting;
  static double baseProbMosqSurvivalResting;
};

#endif
