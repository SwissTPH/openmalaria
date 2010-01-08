/*

  This file is part of OpenMalaria.
 
  Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "WithinHost/WithinHostModel.h"
#include "WithinHost/Descriptive.h"
#include "WithinHost/DescriptiveIPT.h"
#include "WithinHost/Dummy.h"
#include "WithinHost/Empirical.h"
#include "inputData.h"
#include "util/gsl.h"
#include "util/ModelOptions.hpp"
#include "util/errors.hpp"

#include <cmath>
using namespace std;

namespace OM { namespace WithinHost {
    // NOTE: I'd rather use these than x.99 values, but it changes things!
    //const double WithinHostModel::agemin[nages] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 30, 40, 50, 60 };
const double WithinHostModel::agemax[nages] = { 0.99, 1.99, 2.99, 3.99, 4.99, 5.99, 6.99, 7.99, 8.99, 9.99, 10.99, 11.99, 12.99, 13.99, 14.99, 19.99, 24.99, 29.99, 39.99, 49.99, 59.99, 60.99 };
// weight proportions, used by drug code
const double WithinHostModel::wtprop[nages] = { 0.116547265, 0.152531009, 0.181214575, 0.202146126, 0.217216287, 0.237405732, 0.257016899, 0.279053187, 0.293361286, 0.309949502, 0.334474135, 0.350044993, 0.371144279, 0.389814144, 0.412366341, 0.453, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
double WithinHostModel::sigma_i;
double WithinHostModel::immPenalty_22;
double WithinHostModel::asexImmRemain;
double WithinHostModel::immEffectorRemain;
double WithinHostModel::detectionLimit;

// -----  static functions  -----

void WithinHostModel::init() {
  Infection::init();
  sigma_i=sqrt(InputData.getParameter(Params::SIGMA_I_SQ));
  immPenalty_22=1-exp(InputData.getParameter(Params::IMMUNITY_PENALTY));
  immEffectorRemain=exp(-InputData.getParameter(Params::IMMUNE_EFFECTOR_DECAY));
  asexImmRemain=exp(-InputData.getParameter(Params::ASEXUAL_IMMUNITY_DECAY));
  
  double densitybias;
  if (util::ModelOptions::option (util::GARKI_DENSITY_BIAS)) {
      densitybias=InputData.getParameter(Params::DENSITY_BIAS_GARKI);
  } else {
    if ((InputData.get_analysis_no() >= 22) && (InputData.get_analysis_no() <= 30)) {
	cerr << "Warning: these analysis numbers used to mean use Garki density bias. If you do want to use this, specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong." << endl;
    }
    densitybias=InputData.getParameter(Params::DENSITY_BIAS_NON_GARKI);
  }
  detectionLimit=InputData.get_detectionlimit()*densitybias;
  
  if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL)) {
    DummyInfection::init ();
  } else if (util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL)) {
    EmpiricalInfection::initParameters();	// 1-day timestep check
  } else {
    if (util::ModelOptions::option (util::INCLUDES_PK_PD))
	throw util::xml_scenario_error ("INCLUDES_PK_PD is incompatible with the old within-host model");
    DescriptiveInfection::initParameters ();	// 5-day timestep check
    DescriptiveIPTWithinHost::initParameters();
  }
}

size_t WithinHostModel::getAgeGroup (double age) {
    for (size_t i = 0; i < nages; ++i) {
	if (agemax[i] > age)
	    return i;
    }
    return nages-1;	// final category
}

void WithinHostModel::clear() {
  DescriptiveIPTWithinHost::clearParameters();
  DescriptiveInfection::clearParameters();
}

WithinHostModel* WithinHostModel::createWithinHostModel () {
  if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL)) {
    return new DummyWithinHostModel();
  } else if (util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL)) {
    return new EmpiricalWithinHostModel();
  } else {
    if (DescriptiveIPTWithinHost::iptActive)
      return new DescriptiveIPTWithinHost();
    else
      return new DescriptiveWithinHostModel();
  }
}


// -----  Non-static  -----

WithinHostModel::WithinHostModel () :
    _cumulativeh(0.0), _cumulativeY(0.0), _cumulativeYlag(0.0),
    _MOI(0), totalDensity(0.0), timeStepMaxDensity(0.0)
{
    _innateImmSurvFact = exp(-gsl::rngGauss(0, sigma_i));
}


void WithinHostModel::clearInfections (bool) {
  clearAllInfections();
}

void WithinHostModel::IPTiTreatment (SurveyAgeGroup ageGroup) {
  throw util::xml_scenario_error (string ("Timed IPT treatment when no IPT description is present in interventions"));
}


// -----  immunity  -----

void WithinHostModel::updateImmuneStatus(){
  if (immEffectorRemain < 1){
    _cumulativeh*=immEffectorRemain;
    _cumulativeY*=immEffectorRemain;
  }
  if (asexImmRemain < 1){
    _cumulativeh*=asexImmRemain/
        (1+(_cumulativeh*(1-asexImmRemain)/Infection::cumulativeHstar));
    _cumulativeY*=asexImmRemain/
        (1+(_cumulativeY*(1-asexImmRemain)/Infection::cumulativeYstar));
  }
  _cumulativeYlag = _cumulativeY;
}

void WithinHostModel::immunityPenalisation() {
  _cumulativeY = _cumulativeYlag - immPenalty_22*(_cumulativeY-_cumulativeYlag);
  if (_cumulativeY < 0) {
    _cumulativeY=0.0;
  }
}


// -----  Summarize  -----

void WithinHostModel::summarize (Survey& survey, SurveyAgeGroup ageGroup) {
  int patentInfections;
  int numInfections = countInfections (patentInfections);
  if (numInfections) {
    survey.reportInfectedHosts(ageGroup,1);
    survey.addToInfections(ageGroup, numInfections);
    survey.addToPatentInfections(ageGroup, patentInfections);
  }
  // Treatments in the old ImmediateOutcomes clinical model clear infections immediately
  // (and are applied after calculateDensities()); here we report the last calculated density.
  if (parasiteDensityDetectible()) {
    survey.reportPatentHosts(ageGroup, 1);
    survey.addToLogDensity(ageGroup, log(totalDensity));
  }
}


void WithinHostModel::checkpoint (istream& stream) {
    _innateImmSurvFact & stream;
    _cumulativeh & stream;
    _cumulativeY & stream;
    _cumulativeYlag & stream;
    _MOI & stream;
    totalDensity & stream;
    timeStepMaxDensity & stream;
    
    if (_MOI < 0 || _MOI > MAX_INFECTIONS)
	throw util::checkpoint_error ("_MOI");
}
void WithinHostModel::checkpoint (ostream& stream) {
    _innateImmSurvFact & stream;
    _cumulativeh & stream;
    _cumulativeY & stream;
    _cumulativeYlag & stream;
    _MOI & stream;
    totalDensity & stream;
    timeStepMaxDensity & stream;
}

} }