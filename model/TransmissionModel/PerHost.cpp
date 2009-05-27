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
#include "TransmissionModel/PerHost.h"
#include "GSLWrapper.h"
#include "inputData.h"
#include "TransmissionModel/Vector.h"
#include "summary.h"
#include "intervention.h"

double PerHostTransmission::BaselineAvailabilityShapeParam;
const double PerHostTransmission::susceptibility= 0.702;
double PerHostTransmission::gamma_p; 
double PerHostTransmission::Sinf; 
double PerHostTransmission::Simm; 
double PerHostTransmission::Xstar_p; 
double PerHostTransmission::Estar; 
double PerHostTransmission::InfectionrateShapeParam;

void PerHostTransmission::initParameters () {
  EntoInterventionITN::initParameters();
  EntoInterventionIRS::initParameters();
  
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
    InfectionrateShapeParam = (PerHostTransmission::BaselineAvailabilityShapeParam+1.0) / (r_square_Gamma*PerHostTransmission::BaselineAvailabilityShapeParam - 1.0);
    InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  }
  else if (Global::modelVersion &
    (LOGNORMAL_MASS_ACTION | LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM)) {
    InfectionrateShapeParam = sqrt(r_square_LogNormal - 1.86*pow(PerHostTransmission::BaselineAvailabilityShapeParam, 2));
  InfectionrateShapeParam=std::max(InfectionrateShapeParam, 0.0);
  }
}

PerHostTransmission::PerHostTransmission (TransmissionModel& tm) :
  _cumulativeEIRa(0.0), _pinfected(0.0)
{
  VectorTransmission* vTM = dynamic_cast<VectorTransmission*> (&tm);
  if (vTM) {
    species.resize (vTM->numSpecies);
    for (size_t i = 0; i < vTM->numSpecies; ++i)
      species[i].initialise (vTM->species[i]);
  }
  
  //TODO: as for static counterpart
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    _BaselineAvailabilityToMosquitoes=(W_GAMMA((BaselineAvailabilityShapeParam), (BaselineAvailabilityMean/BaselineAvailabilityShapeParam)));
  }
  else if(Global::modelVersion & LOGNORMAL_MASS_ACTION) {
    _BaselineAvailabilityToMosquitoes=(W_LOGNORMAL((log(BaselineAvailabilityMean))-(0.5*pow(BaselineAvailabilityShapeParam, 2)), (BaselineAvailabilityShapeParam)));
  }
  else if (Global::modelVersion & TRANS_HET) {
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

PerHostTransmission::PerHostTransmission (istream& in, TransmissionModel& tm) {
  in >> _cumulativeEIRa;
  in >> _pinfected;
  in >> _BaselineAvailabilityToMosquitoes;
  VectorTransmission* vTM = dynamic_cast<VectorTransmission*> (&tm);
  if (vTM) {
    species.resize (vTM->numSpecies);
    for (vector<HostMosquitoInteraction>::iterator hMI = species.begin(); hMI != species.end(); ++hMI)
      hMI->read (in);
  }
}

void PerHostTransmission::write (ostream& out) const {
  out << _cumulativeEIRa << endl; 
  out << _pinfected << endl; 
  out << _BaselineAvailabilityToMosquitoes << endl; 
  for (vector<HostMosquitoInteraction>::const_iterator hMI = species.begin(); hMI != species.end(); ++hMI)
    hMI->write (out);
}


void PerHostTransmission::summarize (Summary& summary, double age) {
  summary.addToExpectedInfected(age, _pinfected);
}


double PerHostTransmission::getExpectedNumberOfInfections (Human& human, double age_adj_EIR) {
  double baseAvailToMos = _BaselineAvailabilityToMosquitoes;
  //The age-adjusted EIR, possibly adjusted for bed nets.
  double expectedNumInfections;
  
  double ExpectedInfectionRate = age_adj_EIR * baseAvailToMos * susceptibility * Global::interval;
  if (Global::modelVersion & NEGATIVE_BINOMIAL_MASS_ACTION) {
    expectedNumInfections = (W_GAMMA((InfectionrateShapeParam), (ExpectedInfectionRate/InfectionrateShapeParam)));
  } else if (Global::modelVersion & LOGNORMAL_MASS_ACTION) {
    expectedNumInfections = sampleFromLogNormal(W_UNIFORM(),
        log(ExpectedInfectionRate) - 0.5*pow(InfectionrateShapeParam, 2),
        InfectionrateShapeParam);
  } else {
    //The default model is that in Smith et al, AJTMH 2006 75 Suppl 2
    double survivalOfInoculum=(1.0+pow((_cumulativeEIRa/Xstar_p), gamma_p));
    survivalOfInoculum = Simm+(1.0-Simm)/survivalOfInoculum;
    survivalOfInoculum = survivalOfInoculum*(Sinf+(1-Sinf)/(1 + age_adj_EIR/Estar));
    
    if(Global::modelVersion & LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM) {
      expectedNumInfections = survivalOfInoculum *
          sampleFromLogNormal(W_UNIFORM(),
                              log(ExpectedInfectionRate) - 0.5*pow(InfectionrateShapeParam, 2),
                              InfectionrateShapeParam);
    } else {
      expectedNumInfections = survivalOfInoculum *
          age_adj_EIR * Global::interval * baseAvailToMos;
    }
  }
  
  //Introduce the effect of vaccination. Note that this does not affect cumEIR.
  if (Vaccine::PEV.active) {
    expectedNumInfections *= (1 - human.getPEVEfficacy());
  }
  return expectedNumInfections;
}

int PerHostTransmission::numNewInfections (double expectedInfectionRate, double expectedNumberOfInfections){
  //TODO: this code does not allow for variations in baseline availability
  //this is only likely to be relevant in some models but should not be
  //forgotten
  
  //Update pre-erythrocytic immunity
  if (Global::modelVersion & 
    (TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET)) {
    _cumulativeEIRa+=double(Global::interval)*expectedInfectionRate*_BaselineAvailabilityToMosquitoes;
  }
  else {
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


double PerHostTransmission::entoAvailability (size_t speciesIndex) const {
  return species[speciesIndex].entoAvailability
  * species[speciesIndex].entoInterventionITN.availability()
  * species[speciesIndex].entoInterventionIRS.availability();
}
double PerHostTransmission::probMosqBiting (size_t speciesIndex) const {
  return species[speciesIndex].probMosqBiting
  * species[speciesIndex].entoInterventionITN.probMosqBiting();
}
double PerHostTransmission::probMosqFindRestSite (size_t speciesIndex) const {
  return species[speciesIndex].probMosqFindRestSite
  * species[speciesIndex].entoInterventionITN.probMosqFindRestSite();
}
double PerHostTransmission::probMosqSurvivalResting (size_t speciesIndex) const {
  return species[speciesIndex].probMosqSurvivalResting
  * species[speciesIndex].entoInterventionIRS.probMosqSurvivalResting();
}


void HostMosquitoInteraction::initialise (VectorTransmissionSpecies base)
{
  //FIXME: should be partially random:
  entoAvailability = base.entoAvailability;
  probMosqBiting = base.probMosqBiting;
  probMosqFindRestSite = base.probMosqFindRestSite;
  probMosqSurvivalResting = base.probMosqSurvivalResting;
}

void HostMosquitoInteraction::read (istream& in) {
  in >> entoAvailability;
  in >> probMosqBiting;
  in >> probMosqFindRestSite;
  in >> probMosqSurvivalResting;
  in >> entoInterventionITN;
  in >> entoInterventionIRS;
}

void HostMosquitoInteraction::write (ostream& out) const {
  out << entoAvailability << endl;
  out << probMosqBiting << endl;
  out << probMosqFindRestSite << endl;
  out << probMosqSurvivalResting << endl;
  out << entoInterventionITN;
  out << entoInterventionIRS;
}
