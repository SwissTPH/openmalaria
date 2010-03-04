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
#include "WithinHost/Common.h"
#include "WithinHost/DummyInfection.h"
#include "WithinHost/EmpiricalInfection.h"
#include "inputData.h"
#include "util/random.h"
#include "util/ModelOptions.hpp"
#include "util/errors.hpp"

#include <cmath>
using namespace std;
using namespace OM::util;

namespace OM { namespace WithinHost {
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
    if ((InputData().getAnalysisNo() >= 22) && (InputData().getAnalysisNo() <= 30)) {
	cerr << "Warning: these analysis numbers used to mean use Garki density bias. If you do want to use this, specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong." << endl;
    }
    densitybias=InputData.getParameter(Params::DENSITY_BIAS_NON_GARKI);
  }
  detectionLimit=InputData().getMonitoring().getSurveys().getDetectionLimit()*densitybias;
  
  if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL)) {
    DummyInfection::initParameters ();
  } else if (util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL)) {
    EmpiricalInfection::initParameters();	// 1-day timestep check
  } else {
    DescriptiveInfection::initParameters ();	// 5-day timestep check
    DescriptiveIPTWithinHost::initParameters();
  }
}

void WithinHostModel::clear() {
  DescriptiveIPTWithinHost::clearParameters();
  DescriptiveInfection::clearParameters();
}

WithinHostModel* WithinHostModel::createWithinHostModel () {
  if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL) || util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL)) {
    return new CommonWithinHost();
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
    _innateImmSurvFact = exp(-random::gauss(0, sigma_i));
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