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
#ifndef Hmod_NonVectorTransmission
#define Hmod_NonVectorTransmission

#include "TransmissionModel.h"

namespace scnXml {
  class NonVector;
}

//! Base transmission model, as used in Phase A
class NonVectorTransmission : public TransmissionModel { 
public:

  NonVectorTransmission(const scnXml::NonVector& nonVectorData);
  virtual ~NonVectorTransmission();
  //! initialise the main simulation 
  void initMainSimulation (int populationSize);
  
  /** Change the scnXml::NonVector data (changeEIR intervention). */
  void setNonVectorData (const scnXml::NonVector& nonVectorData);

 /*! 
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
  
  /** @brief Initialization function, setting up EIR arrays
   *
   * Reads in the estimates of the EIR for each village and each day
   * and converts this into EIR estimates per five day period
   * assuming that the annual cycle repeated during the pre-intervention period
   */
  void inputEIR (const scnXml::NonVector& nonVectorData); 

  /** Calculates EIR (in adults) during the main period of the simulation,
   * based on vectorial capacity or looks up EIR in the input data.
   * 
   * \param simulationTime Time since start of simulation . */
  virtual double calculateEIR(int simulationTime, PerHostTransmission&); 
 
private:

  /// Processes each daily EIR estimate, allocating each day in turn to the
  /// appropriate time period. EIRdaily is the value of the daily EIR read in
  /// from the .XML file.
  void updateEIR (int day, double EIRdaily); 
  double averageEIR (const scnXml::NonVector& nonVectorData); 
  
//! multiplier used to calculate a positive EIR value where the measured value is zero
/* 
  0.01 was old pv(30) Now a constant. min_EIR_mult multiplies the average EIR to obtain a value used for the EIR during periods when it is too low 
  to be measureable. The value of 0.01 was old pv(30) Now a constant. 
  0.01 was old pv(30) Now a constant. 
*/ 
  static const double min_EIR_mult; 

//! The maximum number of daily EIR values specified for intervention phase 
 /*! 
  The EIR is input on a daily basis even when the interval length is >1 day 
  Average EIR per interval is computed in initEntoParameters 
  \sa initEntoParameters 
  */ 
  static const int maxDurIntPhaseEIR= 1500; 
  
  /** @brief Variables set by constructor.
   *
   * There shouldn't be any need to checkpoint these, at least before
   * interventions take effect. */
  //@{
  //! The maximum number of intervals in the intervention phase.
  int maxIntervals; 
  
  //! Number of days contributing to each EIR estimate for pre-intervention 
  int *nDays;

  /// Number of days contributing to each EIR estimate (post intervention) 
  int *ino;

  ///intEIR() EIR per time interval during the intervention period 
  double *intEIR;
  
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
  double gamma_p; 
 
//!Lower limit of success probability of inoculations at high exposure in Phase A model 
 /*! 
  \sa getExpectedNumberOfInfections(),gamma_p,Simm,Xstar_p,Estar 
  */ 
  double Sinf; 
 
//!Lower limit of success probability of inoculations in immune individuals in Phase A model 
 /*! 
  \sa getExpectedNumberOfInfections(),gamma_p,Sinf,Xstar_p,Estar 
  */ 
  double Simm; 
 
//!Critical value of cumulative number of entomologic inoculations in Phase A model 
 /*! 
  \sa getExpectedNumberOfInfections(),gamma_p,Sinf,Simm,Estar 
  */ 
  double Xstar_p; 
 
//!Critical value of EIR in Phase A pre-erythrocytic model 
 /*! 
  \sa getExpectedNumberOfInfections(),gamma_p,Sinf,Simm,Xstar_p 
  */ 
  double Estar; 
 

//! Describes the shape of the Infectionrate distribution, related to the baseline availabilty distr. 
  double InfectionrateShapeParam;

//! Variance of Infection Rate according to fielddata 
  static const double totalInfectionrateVariance; 
  
  //! The duration of sporogony in time steps 
  int nspore;
  //@}
};
#endif
