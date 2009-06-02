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

#include "InfectionIncidenceModel.h"
#include "GSLWrapper.h"
#include "inputData.h"
#include "human.h"	//FIXME: remove

double InfectionIncidenceModel::BaselineAvailabilityShapeParam;
const double InfectionIncidenceModel::susceptibility= 0.702;
double InfectionIncidenceModel::gamma_p; 
double InfectionIncidenceModel::Sinf; 
double InfectionIncidenceModel::Simm; 
double InfectionIncidenceModel::Xstar_p; 
double InfectionIncidenceModel::Estar; 
double InfectionIncidenceModel::InfectionrateShapeParam;

// -----  static initialisation  -----

void InfectionIncidenceModel::init () {
  BaselineAvailabilityShapeParam=getParameter(Params::BASELINE_AVAILABILITY_SHAPE);
  
  // NOTE: The following concerns translating an EIR into a number of
  // infections, and could be placed in a different module.
  // TODO: use class inheritance to replace if blocks
  
  gamma_p=getParameter(Params::GAMMA_P);
  Sinf=1-exp(-getParameter(Params::NEG_LOG_ONE_MINUS_SINF));
  Simm=getParameter(Params::SIMM);
  Estar=getParameter(Params::E_STAR);
  Xstar_p=getParameter(Params::X_STAR_P);
  
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
  
  //! constant defining the constraint for the log Normal variance
  /// Used for the case where availability is assumed log Normally distributed
  double r_square_LogNormal = log(1.0+r_square_Gamma);
  
  //TODO: Sanity check for sqrt and division by zero
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    InfectionrateShapeParam = (BaselineAvailabilityShapeParam+1.0) / (r_square_Gamma*BaselineAvailabilityShapeParam - 1.0);
    InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  }
  else if (Global::modelVersion &
    (LOGNORMAL_MASS_ACTION | LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM)) {
    InfectionrateShapeParam = sqrt(r_square_LogNormal - 1.86*pow(BaselineAvailabilityShapeParam, 2));
  InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  }
}


// -----  non-static non-checkpointing constructors  -----

InfectionIncidenceModel* InfectionIncidenceModel::createModel () {
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    return new NegBinomMAII ();
  } else if(Global::modelVersion & LOGNORMAL_MASS_ACTION) {
    return new LogNormalMAII ();
  } else if(Global::modelVersion & LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM) {
    return new LogNormalMAPlusPreImmII();
  } else
    return new InfectionIncidenceModel();
}

InfectionIncidenceModel::InfectionIncidenceModel () :
  _pinfected(0.0), _cumulativeEIRa(0.0)
{
  if (Global::modelVersion & TRANS_HET) {
    _BaselineAvailabilityToMosquitoes=0.2;
    if (W_UNIFORM() < 0.5) {
      _BaselineAvailabilityToMosquitoes=1.8;
    }
  }
  else {
    _BaselineAvailabilityToMosquitoes=BaselineAvailabilityMean;
  }
  // NOTE: _BaselineAvailabilityToMosquitoes MAY be re-set in the Human constructor
}
InfectionIncidenceModel::InfectionIncidenceModel (double bATM) :
  _pinfected(0.0), _BaselineAvailabilityToMosquitoes(bATM), _cumulativeEIRa(0.0)
{
  // NOTE: _BaselineAvailabilityToMosquitoes MAY be re-set in the Human constructor
}
NegBinomMAII::NegBinomMAII () :
  InfectionIncidenceModel(W_GAMMA(BaselineAvailabilityShapeParam,
				  BaselineAvailabilityMean/BaselineAvailabilityShapeParam))
{}
LogNormalMAII::LogNormalMAII () :
  InfectionIncidenceModel(W_LOGNORMAL(log(BaselineAvailabilityMean)-(0.5*pow(BaselineAvailabilityShapeParam, 2)),
				       BaselineAvailabilityShapeParam))
{}
LogNormalMAPlusPreImmII::LogNormalMAPlusPreImmII () :
  InfectionIncidenceModel()
{}


// -----  checkpointing  -----

InfectionIncidenceModel* InfectionIncidenceModel::createModel (istream& in) {
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    return new NegBinomMAII (in);
  } else if(Global::modelVersion & LOGNORMAL_MASS_ACTION) {
    return new LogNormalMAII (in);
  } else if(Global::modelVersion & LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM) {
    return new LogNormalMAPlusPreImmII(in);
  } else
    return new InfectionIncidenceModel(in);
}

InfectionIncidenceModel::InfectionIncidenceModel (istream& in) {
  in >> _cumulativeEIRa;
  in >> _pinfected;
  in >> _BaselineAvailabilityToMosquitoes;
}
NegBinomMAII::NegBinomMAII (istream& in) :
  InfectionIncidenceModel(in)
{}
LogNormalMAII::LogNormalMAII (istream& in) :
  InfectionIncidenceModel(in)
{}
LogNormalMAPlusPreImmII::LogNormalMAPlusPreImmII (istream& in) :
  InfectionIncidenceModel(in)
{}

void InfectionIncidenceModel::write (ostream& out) const {
  out << _cumulativeEIRa << endl; 
  out << _pinfected << endl; 
  out << _BaselineAvailabilityToMosquitoes << endl; 
}


// -----  non-static methods  -----

void InfectionIncidenceModel::summarize (Summary& summary, double age) {
  summary.addToExpectedInfected(age, _pinfected);
}


double InfectionIncidenceModel::getExpectedNumberOfInfections (Human& human, double age_adj_EIR) {
  double expectedNumInfections = getModelExpectedInfections (age_adj_EIR);
  
  //Introduce the effect of vaccination. Note that this does not affect cumEIR.
  if (Vaccine::PEV.active) {
    expectedNumInfections *= (1 - human.getPEVEfficacy());
  }
  return expectedNumInfections;
}

double InfectionIncidenceModel::getModelExpectedInfections (double age_adj_EIR) {
  return survivalOfInoculum(age_adj_EIR) * age_adj_EIR * Global::interval * _BaselineAvailabilityToMosquitoes;
}
double NegBinomMAII::getModelExpectedInfections (double age_adj_EIR) {
  double ExpectedInfectionRate = age_adj_EIR * _BaselineAvailabilityToMosquitoes * susceptibility * Global::interval;
  return (W_GAMMA((InfectionrateShapeParam), (ExpectedInfectionRate/InfectionrateShapeParam)));
}
double LogNormalMAII::getModelExpectedInfections (double age_adj_EIR) {
  double ExpectedInfectionRate = age_adj_EIR * _BaselineAvailabilityToMosquitoes * susceptibility * Global::interval;
  return sampleFromLogNormal(W_UNIFORM(),
			     log(ExpectedInfectionRate) - 0.5*pow(InfectionrateShapeParam, 2),
			     InfectionrateShapeParam);
}
double LogNormalMAPlusPreImmII::getModelExpectedInfections (double age_adj_EIR) {
  double ExpectedInfectionRate = age_adj_EIR * _BaselineAvailabilityToMosquitoes * susceptibility * Global::interval;
  
  return survivalOfInoculum(age_adj_EIR) *
    sampleFromLogNormal(W_UNIFORM(),
			log(ExpectedInfectionRate) - 0.5*pow(InfectionrateShapeParam, 2),
			InfectionrateShapeParam);
}

double InfectionIncidenceModel::survivalOfInoculum (double age_adj_EIR) {
  double survivalOfInoculum=(1.0+pow((_cumulativeEIRa/Xstar_p), gamma_p));
  survivalOfInoculum = Simm+(1.0-Simm)/survivalOfInoculum;
  survivalOfInoculum = survivalOfInoculum*(Sinf+(1-Sinf)/(1 + age_adj_EIR/Estar));
  return survivalOfInoculum;
}

int InfectionIncidenceModel::numNewInfections (double expectedInfectionRate, double expectedNumberOfInfections){
  //TODO: this code does not allow for variations in baseline availability
  //this is only likely to be relevant in some models but should not be
  //forgotten
  
  //Update pre-erythrocytic immunity
  if (Global::modelVersion & 
    (TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET)) {
    _cumulativeEIRa+=double(Global::interval)*expectedInfectionRate*_BaselineAvailabilityToMosquitoes;
  } else {
    _cumulativeEIRa+=double(Global::interval)*expectedInfectionRate;
  }
  
  _pinfected = 1.0 - exp(-expectedNumberOfInfections) * (1.0-_pinfected);
  if (_pinfected < 0.0)
    _pinfected = 0.0;
  else if (_pinfected > 1.0)
    _pinfected = 1.0;
  
  if (expectedNumberOfInfections > 0.0000001)
    return W_POISSON(expectedNumberOfInfections);
  else
    return 0;
}
