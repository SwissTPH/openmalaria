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

#ifndef Hmod_InfectionIncidence
#define Hmod_InfectionIncidence

#include "global.h"
#include "summary.h"
#include "inputData.h"
#include "intervention.h"

class Human;	//FIXME: remove this

/** Models how a per-host EIR translates into new infections
 * (roughly when bites from infected mosquitos infect the host).
 *
 * There are four versions of this model:
 * - InfectionIncidenceModel (the default): Smith et al, AJTMH 2006 75 Suppl 2
 * - NegBinomMAII: NEGATIVE_BINOMIAL_MASS_ACTION
 * - LogNormalMAII: LOGNORMAL_MASS_ACTION
 * - LogNormalMAPlusPreImmII: LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM
 */
class InfectionIncidenceModel
{
public:
  ///@brief Static initialisation & constructors
  //@{
  /// Read in/initialise parameters
  static void init();
  /// Create a new instance
  static InfectionIncidenceModel* createModel ();
  /// Read from a checkpoint
  static InfectionIncidenceModel* createModel (istream& in);
  //@}
  
protected:
  /// Create a new model, deciding what _BaselineAvailabilityToMosquitoes to use
  InfectionIncidenceModel ();
  /// Create a new model, passing the _BaselineAvailabilityToMosquitoes from
  /// a derived class's constructor
  InfectionIncidenceModel (double bATM);
  /// Load from a checkpoint
  InfectionIncidenceModel (istream& in);
  
public:
  /// Write a checkpoint
  void write (ostream& out) const;
  
  /// Output _pinfected to the summary
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
  double getExpectedNumberOfInfections (Human& human, double age_adj_EIR);
  
  //! Calculate the number of new infections to introduce via a stochastic process
  int numNewInfections(double expectedInfectionRate, double expectedNumberOfInfections);
  
protected:
  /// Calculates the expected number of infections, excluding vaccine effects
  virtual double getModelExpectedInfections (double age_adj_EIR);
  
  double survivalOfInoculum (double age_adj_EIR);
  
public:	//TODO - maybe better if not public
  /** Probability of infection (cumulative or reset to zero in massTreatment).
   *
   * Appears to be used only for calculating expected inoculations for the
   * analysis of pre-erythrocytic immunity. */
  double _pinfected;
  /// Baseline availability to mosquitoes
  double _BaselineAvailabilityToMosquitoes;
private:
  //!Number of infective bites since birth
    double _cumulativeEIRa;//TODO: not needed by NegBinomMAII and LogNormalMAII
  
protected:	// Static data
  /* Shape constant of (Gamma) distribution of availability
  real, parameter :: BaselineAvailabilityGammaShapeParam =1.0 */
  static double BaselineAvailabilityShapeParam;
  
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

class NegBinomMAII : public InfectionIncidenceModel {
public:
  NegBinomMAII ();
  NegBinomMAII (istream& in);
protected:
  double getModelExpectedInfections (double age_adj_EIR);
};
class LogNormalMAII : public InfectionIncidenceModel {
public:
  LogNormalMAII ();
  LogNormalMAII (istream& in);
protected:
  double getModelExpectedInfections (double age_adj_EIR);
};
class LogNormalMAPlusPreImmII : public InfectionIncidenceModel {
public:
  LogNormalMAPlusPreImmII ();
  LogNormalMAPlusPreImmII (istream& in);
protected:
  double getModelExpectedInfections (double age_adj_EIR);
};

#endif