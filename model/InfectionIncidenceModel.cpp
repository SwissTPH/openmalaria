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

#include "global.h"
#include "InfectionIncidenceModel.h"
#include "GSLWrapper.h"
#include "inputData.h"
#include "Transmission/PerHost.h"

double InfectionIncidenceModel::BaselineAvailabilityShapeParam;
double InfectionIncidenceModel::InfectionrateShapeParam;

double InfectionIncidenceModel::gamma_p;
double InfectionIncidenceModel::Sinf;
double InfectionIncidenceModel::Simm;
double InfectionIncidenceModel::Xstar_pInv;
double InfectionIncidenceModel::EstarInv;

// -----  static initialisation  -----

void InfectionIncidenceModel::init () {
  BaselineAvailabilityShapeParam=getParameter(Params::BASELINE_AVAILABILITY_SHAPE);
  
  gamma_p=getParameter(Params::GAMMA_P);
  Sinf=1-exp(-getParameter(Params::NEG_LOG_ONE_MINUS_SINF));
  Simm=getParameter(Params::SIMM);
  EstarInv = 1.0/getParameter(Params::E_STAR);
  Xstar_pInv = 1.0/getParameter(Params::X_STAR_P);
  
  //TODO: Sanity check for sqrt and division by zero
  //! constant defining the constraint for the Gamma shape parameters
  /// Used for the case where availability is assumed gamma distributed
  double r_square_Gamma;
  /*
  //! Expected number of inoculations
  /// Product of measured EIR, susceptibility and length of time Global::interval
  double gsi = 1.0;
  
  r_square_Gamma=(totalInfectionrateVariance**2-gsi*BaselineAvailabilityMean)/(gsi*BaselineAvailabilityMean)**2
  r_square_Gamma must be greater than zero, so r_square_LogNormal is also. 
  */
  r_square_Gamma=0.649;
  //such that r_square_LogNormal =0.5
  
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    InfectionrateShapeParam = (BaselineAvailabilityShapeParam+1.0) / (r_square_Gamma*BaselineAvailabilityShapeParam - 1.0);
    InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  } else if (Global::modelVersion & LOGNORMAL_MASS_ACTION) {
    //! constant defining the constraint for the log Normal variance
    /// Used for the case where availability is assumed log Normally distributed
    double r_square_LogNormal = log(1.0+r_square_Gamma);
  
    InfectionrateShapeParam = sqrt(r_square_LogNormal - 1.86*pow(BaselineAvailabilityShapeParam, 2));
    InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  }
  
  if (Global::modelVersion & ANY_TRANS_HET)
    cerr << "Warning: will use heterogeneity workaround." << endl;
}


// -----  non-static non-checkpointing constructors  -----

InfectionIncidenceModel* InfectionIncidenceModel::createModel () {
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    return new NegBinomMAII ();
  } else if(Global::modelVersion & LOGNORMAL_MASS_ACTION) {
    return new LogNormalMAII ();
  } else {
    if (Global::modelVersion & ANY_TRANS_HET)
      return new HeterogeneityWorkaroundII();
    else
      return new InfectionIncidenceModel();
  }
}

InfectionIncidenceModel::InfectionIncidenceModel () :
  _pinfected(0.0), _cumulativeEIRa(0.0)
{}


// -----  checkpointing  -----

InfectionIncidenceModel* InfectionIncidenceModel::createModel (istream& in) {
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    return new NegBinomMAII (in);
  } else if(Global::modelVersion & LOGNORMAL_MASS_ACTION) {
    return new LogNormalMAII (in);
  } else {
    if (Global::modelVersion & ANY_TRANS_HET)
      return new HeterogeneityWorkaroundII(in);
    else
      return new InfectionIncidenceModel(in);
  }
}

InfectionIncidenceModel::InfectionIncidenceModel (istream& in) {
  in >> _cumulativeEIRa;
  in >> _pinfected;
}
NegBinomMAII::NegBinomMAII (istream& in) :
  InfectionIncidenceModel(in)
{}
LogNormalMAII::LogNormalMAII (istream& in) :
  InfectionIncidenceModel(in)
{}

void InfectionIncidenceModel::write (ostream& out) const {
  out << _cumulativeEIRa << endl; 
  out << _pinfected << endl; 
}


// -----  non-static methods  -----

double InfectionIncidenceModel::getAvailabilityFactor(double baseAvailability) {
  return baseAvailability;
}
double NegBinomMAII::getAvailabilityFactor(double baseAvailability) {
  return W_GAMMA(BaselineAvailabilityShapeParam,
		 baseAvailability/BaselineAvailabilityShapeParam);
}
double LogNormalMAII::getAvailabilityFactor(double baseAvailability) {
  return W_LOGNORMAL(log(baseAvailability)-(0.5*pow(BaselineAvailabilityShapeParam, 2)),
		     BaselineAvailabilityShapeParam);
}

void InfectionIncidenceModel::summarize (Summary& summary, double age) {
  summary.addToExpectedInfected(age, _pinfected);
}


double InfectionIncidenceModel::getModelExpectedInfections (double effectiveEIR, PerHostTransmission&) {
  // First two lines are availability adjustment: S_1(i,t) from AJTMH 75 (suppl 2) p12 eqn. (5)
  return (Sinf+(1-Sinf) / 
    (1 + effectiveEIR/Global::interval*EstarInv)) *
    susceptibility() * effectiveEIR;
}
double HeterogeneityWorkaroundII::getModelExpectedInfections (double effectiveEIR, PerHostTransmission& phTrans) {
  return (Sinf+(1-Sinf)/(1 + effectiveEIR/(Global::interval*phTrans.entoAvailability())*EstarInv)) *
    susceptibility() * effectiveEIR;
}
double NegBinomMAII::getModelExpectedInfections (double effectiveEIR, PerHostTransmission&) {
  return W_GAMMA(InfectionrateShapeParam,
		  effectiveEIR * susceptibility() / InfectionrateShapeParam);
}
double LogNormalMAII::getModelExpectedInfections (double effectiveEIR, PerHostTransmission&) {
  return sampleFromLogNormal(W_UNIFORM(),
			     log(effectiveEIR * susceptibility()) - 0.5*pow(InfectionrateShapeParam, 2),
			     InfectionrateShapeParam);
}

double InfectionIncidenceModel::susceptibility () {
  if (Global::modelVersion & NO_PRE_ERYTHROCYTIC) {
    //! The average proportion of bites from sporozoite positive mosquitoes resulting in infection. 
    /*! 
    This is computed as 0.19 (the value S from a neg bin mass action model fitted 
    to Saradidi data, divided by 0.302 (the ratio of body surface area in a 
    0.5-6 year old child (as per Saradidi) to adult) 
    \sa getExpectedNumberOfInfections() 
    */ 
    return 0.702;
  } else {
    // S_2(i,t) from AJTMH 75 (suppl 2) p12 eqn. (7)
    return Simm + (1.0-Simm) /
      (1.0 + pow(_cumulativeEIRa*Xstar_pInv, gamma_p));
  }
}

int InfectionIncidenceModel::numNewInfections (double effectiveEIR, double PEVEfficacy, PerHostTransmission& phTrans) {
  double expectedNumInfections = getModelExpectedInfections (effectiveEIR, phTrans);
  //NOTE: error check (should be OK if kappa is checked, for nonVector model)
  if (!finite(effectiveEIR)) {
    ostringstream out;
    out << "Error: effectiveEIR is not finite: " << effectiveEIR << endl;
    throw overflow_error (out.str());
  }
  
  //Introduce the effect of vaccination. Note that this does not affect cumEIR.
  if (Vaccine::PEV.active) {
    expectedNumInfections *= (1.0 - PEVEfficacy);
  }
  
  //Update pre-erythrocytic immunity
  _cumulativeEIRa+=effectiveEIR;
  
  _pinfected = 1.0 - exp(-expectedNumInfections) * (1.0-_pinfected);
  if (_pinfected < 0.0)
    _pinfected = 0.0;
  else if (_pinfected > 1.0)
    _pinfected = 1.0;
  
  if (expectedNumInfections > 0.0000001)
    return W_POISSON(expectedNumInfections);
  else
    return 0;
}
