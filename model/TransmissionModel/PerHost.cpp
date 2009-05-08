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

#include "summary.h"

double PerHostTransmission::BaselineAvailabilityShapeParam;
double PerHostTransmission::baseEntoAvailability;
double PerHostTransmission::baseProbMosqSurvivalBiting;
double PerHostTransmission::baseProbMosqSurvivalResting;

void PerHostTransmission::initParameters () {
  EntoInterventionITN::initParameters();
  EntoInterventionIRS::initParameters();
  
  BaselineAvailabilityShapeParam=getParameter(Params::BASELINE_AVAILABILITY_SHAPE);
  
  baseEntoAvailability = 1.0;		// FIXME: get from xml
  baseProbMosqSurvivalBiting = 1.0;	// FIXME: get from xml
  baseProbMosqSurvivalResting = 1.0;	// FIXME: get from xml
}

PerHostTransmission::PerHostTransmission () :
  _cumulativeEIRa(0.0), _pinfected(0.0)
{
  //FIXME: should be partially random:
  _entoAvailability = baseEntoAvailability;
  _probMosqSurvivalBiting = baseProbMosqSurvivalBiting;
  _probMosqSurvivalResting = baseProbMosqSurvivalResting;
  
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

void PerHostTransmission::read (istream& in) {
  in >> _cumulativeEIRa; 
  in >> _pinfected; 
  in >> _BaselineAvailabilityToMosquitoes; 
  in >> _entoAvailability;
  in >> _probMosqSurvivalBiting;
  in >> _probMosqSurvivalResting;
  in >> entoInterventionITN;
  in >> entoInterventionIRS;
}

void PerHostTransmission::write (ostream& out) const {
  out << _cumulativeEIRa << endl; 
  out << _pinfected << endl; 
  out << _BaselineAvailabilityToMosquitoes << endl; 
  out << _entoAvailability << endl;
  out << _probMosqSurvivalBiting << endl;
  out << _probMosqSurvivalResting << endl;
  out << entoInterventionITN;
  out << entoInterventionIRS;
}


void PerHostTransmission::summarize (Summary& summary, double age) {
  summary.addToExpectedInfected(age, _pinfected);
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

double PerHostTransmission::entoAvailability () const {
  return _entoAvailability
  * entoInterventionITN.availability()
  * entoInterventionIRS.availability();
}
double PerHostTransmission::probMosqSurvivalBiting () const {
  return _probMosqSurvivalBiting
  * entoInterventionITN.probMosqBiting();
}
double PerHostTransmission::probMosqSurvivalResting () const {
  return _probMosqSurvivalResting
  * entoInterventionIRS.probMosqSurvivalResting();
}
