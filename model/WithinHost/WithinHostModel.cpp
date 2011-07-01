/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "WithinHost/WithinHostModel.h"
#include "WithinHost/DescriptiveWithinHost.h"
#include "WithinHost/DescriptiveIPTWithinHost.h"
#include "WithinHost/CommonWithinHost.h"
#include "WithinHost/Infection/DummyInfection.h"
#include "WithinHost/Infection/EmpiricalInfection.h"
#include "WithinHost/Infection/MolineauxInfection.h"
#include "WithinHost/Infection/PennyInfection.h"
#include "inputData.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
//using namespace std;

#include <cmath>
#include <boost/format.hpp>


namespace OM { namespace WithinHost {

using namespace OM::util;

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
    DummyInfection::init ();
  } else if (util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL)) {
    EmpiricalInfection::init();	// 1-day timestep check
  } else if (util::ModelOptions::option (util::MOLINEAUX_WITHIN_HOST_MODEL)) {
    MolineauxInfection::init();
  } else if (util::ModelOptions::option (util::PENNY_WITHIN_HOST_MODEL)) {
      PennyInfection::init();
  } else {
    DescriptiveInfection::init ();	// 5-day timestep check
  }
}

WithinHostModel* WithinHostModel::createWithinHostModel () {
  if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL) ||
      util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL) ||
      util::ModelOptions::option (util::MOLINEAUX_WITHIN_HOST_MODEL) ||
      util::ModelOptions::option (util::PENNY_WITHIN_HOST_MODEL)) {
    return new CommonWithinHost();
  } else {
    if ( util::ModelOptions::option( IPTI_SP_MODEL ) )
      return new DescriptiveIPTWithinHost();
    else
      return new DescriptiveWithinHostModel();
  }
}


// -----  Non-static  -----

WithinHostModel::WithinHostModel () :
    _cumulativeh(0.0), _cumulativeY(0.0), _cumulativeYlag(0.0),
    numInfs(0), totalDensity(0.0), timeStepMaxDensity(0.0)
{
    _innateImmSurvFact = exp(-random::gauss(0, sigma_i));
}


void WithinHostModel::clearInfections (TimeStep now, bool) {
  clearAllInfections();
}

void WithinHostModel::continuousIPT (Monitoring::AgeGroup, bool) {
  throw util::xml_scenario_error (string ("Continuous IPT treatment when no IPT description is present in interventions"));
}
void WithinHostModel::timedIPT (Monitoring::AgeGroup, bool) {
  throw util::xml_scenario_error (string ("Timed IPT treatment when no IPT description is present in interventions"));
}
bool WithinHostModel::hasIPTiProtection (TimeStep maxInterventionAge) const {
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

bool WithinHostModel::summarize (Monitoring::Survey& survey, Monitoring::AgeGroup ageGroup) {
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
    return true;
  }
  return false;
}


void WithinHostModel::checkpoint (istream& stream) {
    _innateImmSurvFact & stream;
    _cumulativeh & stream;
    _cumulativeY & stream;
    _cumulativeYlag & stream;
    numInfs & stream;
    totalDensity & stream;
    timeStepMaxDensity & stream;
    
    if (numInfs > MAX_INFECTIONS)
	throw util::checkpoint_error( (boost::format("numInfs: %1%") %numInfs).str() );
}
void WithinHostModel::checkpoint (ostream& stream) {
    _innateImmSurvFact & stream;
    _cumulativeh & stream;
    _cumulativeY & stream;
    _cumulativeYlag & stream;
    numInfs & stream;
    totalDensity & stream;
    timeStepMaxDensity & stream;
}

} }
