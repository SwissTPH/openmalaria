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

class PerHostTransmission;

/** Models how a per-host EIR translates into new infections
 * (roughly when bites from infected mosquitos infect the host).
 *
 * There are four versions of this model, with different availability models:
 * - InfectionIncidenceModel (the default): Smith et al, AJTMH 2006 75 Suppl 2
 * - HeterogeneityWorkaroundII: emulates old, presumably unintended, behaviour
 * - NegBinomMAII: NEGATIVE_BINOMIAL_MASS_ACTION
 * - LogNormalMAII: LOGNORMAL_MASS_ACTION
 * 
 * There are also two susceptibility models which should be compatible with all
 * of these (see susceptibility()). */
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
  /// Create a new model
  InfectionIncidenceModel ();
  /// Load from a checkpoint
  InfectionIncidenceModel (istream& in);
public:
  virtual ~InfectionIncidenceModel() {}
  
  /// Write a checkpoint
  void write (ostream& out) const;
  
  /** Return an availability multiplier, dependant on the model (NegBinomMAII
   * and LogNormalMAII models use this). Ideally availability adjustments
   * should have nothing to do with the InfectionIncidenceModel though.
   * 
   * @param baseAvailability This was BaselineAvailabilityMean from Constant.h,
   *	and had the value 1.0. Whether it should be anything else I don't know.
   */
  virtual double getAvailabilityFactor(double baseAvailability = 1.0);
  
  /// Output _pinfected to the summary
  void summarize (Summary&, double age);
  
  /** Calculate the number of new infections to introduce.
   * 
   * Firstly converts the EIR into an expected number of infections (what
   * getExpectedNumberOfInfections() used to do):
   * 1. Calculates h from the EIR measured on adults where h is the expected
   * number of epidemiological inoculations.
   * 2. Calculates the updated values of the pre-erythrocytic exposure.
   * 
   * Secondly calculates the number of new infections to introduce via a
   * stochastic process. */
  int numNewInfections(double effectiveEIR, double PEVEfficacy, PerHostTransmission& phTrans);
  
protected:
  /// Calculates the expected number of infections, excluding vaccine effects
  virtual double getModelExpectedInfections (double effectiveEIR, PerHostTransmission& phTrans);
  
  double susceptibility ();
  
  /** Probability of infection (cumulative or reset to zero in massTreatment).
   *
   * Appears to be used only for calculating expected inoculations for the
   * analysis of pre-erythrocytic immunity. */
  double _pinfected;
  
  //!Number of infective bites since birth
  double _cumulativeEIRa;//TODO: not needed by NegBinomMAII and LogNormalMAII
  
protected:	// Static data
  /* Shape constant of (Gamma) distribution of availability
  real, parameter :: BaselineAvailabilityGammaShapeParam =1.0 */
  static double BaselineAvailabilityShapeParam;
  
  //VARIABLES INCLUDED IN CORE GETs of number of infections 
  //! Describes the shape of the Infectionrate distribution, related to the baseline availabilty distr. 
  static double InfectionrateShapeParam;
  
  /** @brief Variables for calculating survivalOfInoculum() */
  //@{
  //!Steepness of relationship between success of inoculation and Xp in Phase A model 
  static double gamma_p;
  //!Lower limit of success probability of inoculations at high exposure in Phase A model 
  static double Sinf;
  //!Lower limit of success probability of inoculations in immune individuals in Phase A model 
  static double Simm;
  //!1 over the critical value of cumulative number of entomologic inoculations in Phase A model 
  static double Xstar_pInv;
  //!1 over the critical value of EIR in Phase A pre-erythrocytic model 
  static double EstarInv;
  //@}
};

/** A workaround to produce the same results as with heterogeneity work-units.
 *
 * The EIR passed into the function was not in one place adjusted by the
 * availability factor used in transmission heterogeneity, where it possibly
 * should have been. In any case, this should allow reproducing those results.
 */
class HeterogeneityWorkaroundII : public InfectionIncidenceModel {
public:
  HeterogeneityWorkaroundII () {}
  HeterogeneityWorkaroundII (istream& in) :
    InfectionIncidenceModel (in) {}
  virtual ~HeterogeneityWorkaroundII() {}
protected:
  double getModelExpectedInfections (double effectiveEIR, PerHostTransmission& phTrans);
};
class NegBinomMAII : public InfectionIncidenceModel {
public:
  NegBinomMAII () {}
  NegBinomMAII (istream& in);
  virtual ~NegBinomMAII() {}
  virtual double getAvailabilityFactor(double baseAvailability = 1.0);
protected:
  double getModelExpectedInfections (double effectiveEIR, PerHostTransmission&);
};
class LogNormalMAII : public InfectionIncidenceModel {
public:
  LogNormalMAII () {}
  LogNormalMAII (istream& in);
  virtual ~LogNormalMAII() {}
  virtual double getAvailabilityFactor(double baseAvailability = 1.0);
protected:
  double getModelExpectedInfections (double effectiveEIR, PerHostTransmission&);
};

#endif
