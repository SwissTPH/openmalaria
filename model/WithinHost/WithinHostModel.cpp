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
#include <stdexcept>

using namespace std;

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

// -----  Initialization  -----

void WithinHostModel::init() {
  sigma_i=sqrt(getParameter(Params::SIGMA_I_SQ));
  immPenalty_22=1-exp(getParameter(Params::IMMUNITY_PENALTY));
  immEffectorRemain=exp(-getParameter(Params::IMMUNE_EFFECTOR_DECAY));
  asexImmRemain=exp(-getParameter(Params::ASEXUAL_IMMUNITY_DECAY));
  
  double densitybias;
  /*
  TODO: This densitiybias function should be part of the scenario description XML, not the parameter element.
  or maybe it should be a parameter, as we want to fit it... but the garki analysis numbers are a bit dangerous
  add an attribute to scenario.xml densityQuantification="malariaTherapy|garki|other"
  */
  if (( get_analysis_no() <  22) || ( get_analysis_no() >  30)) {
    densitybias=getParameter(Params::DENSITY_BIAS_NON_GARKI);
  }
  else {
    densitybias=getParameter(Params::DENSITY_BIAS_GARKI);
  }
  detectionLimit=get_detectionlimit()*densitybias;
  
  if (Global::modelVersion & DUMMY_WITHIN_HOST_MODEL) {
    DummyInfection::init ();
  } else if (Global::modelVersion & EMPIRICAL_WITHIN_HOST_MODEL) {
    EmpiricalInfection::initParameters();	// 1-day timestep check
  } else {
    if (Global::modelVersion & INCLUDES_PK_PD)
      throw xml_scenario_error ("INCLUDES_PK_PD is incompatible with the old within-host model");
    DescriptiveInfection::initParameters ();	// 5-day timestep check
    DescriptiveIPTWithinHost::initParameters();
  }
}

void WithinHostModel::clear() {
  DescriptiveIPTWithinHost::clearParameters();
  DescriptiveInfection::clearParameters();
}

WithinHostModel* WithinHostModel::createWithinHostModel () {
  if (Global::modelVersion & DUMMY_WITHIN_HOST_MODEL) {
    return new DummyWithinHostModel();
  } else if (Global::modelVersion & EMPIRICAL_WITHIN_HOST_MODEL) {
    return new EmpiricalWithinHostModel();
  } else {
    if (DescriptiveIPTWithinHost::iptActive)
      return new DescriptiveIPTWithinHost();
    else
      return new DescriptiveWithinHostModel();
  }
}

WithinHostModel* WithinHostModel::createWithinHostModel (istream& in) {
  if (Global::modelVersion & DUMMY_WITHIN_HOST_MODEL) {
    return new DummyWithinHostModel(in);
  } else if (Global::modelVersion & EMPIRICAL_WITHIN_HOST_MODEL) {
    return new EmpiricalWithinHostModel(in);
  } else {
    if (DescriptiveIPTWithinHost::iptActive)
      return new DescriptiveIPTWithinHost(in);
    else
      return new DescriptiveWithinHostModel(in);
  }
}

WithinHostModel::WithinHostModel(istream& in) {
  in >> totalDensity;
  in >> timeStepMaxDensity;
  in >> _cumulativeh;
  in >> _cumulativeY;
  in >> _cumulativeYlag;
}
void WithinHostModel::write (ostream& out) const {
  out << totalDensity << endl;
  out << timeStepMaxDensity << endl;
  out << _cumulativeh << endl;
  out << _cumulativeY << endl;
  out << _cumulativeYlag << endl;
}

void WithinHostModel::clearInfections (bool) {
  clearAllInfections();
}

void WithinHostModel::IPTiTreatment (int ageGroup) {
  throw xml_scenario_error (string ("Timed IPT treatment when no IPT description is present in interventions"));
}

size_t WithinHostModel::getAgeGroup (double age) {
  for (size_t i = 0; i < nages; ++i) {
    if (agemax[i] > age)
      return i;
  }
  return nages-1;	// final category
}


// -----  immunity  -----

double WithinHostModel::immunitySurvivalFactor () {
}

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
  _cumulativeY=(double)_cumulativeYlag-(immPenalty_22*(_cumulativeY-_cumulativeYlag));
  if (_cumulativeY <  0) {
    _cumulativeY=0.0;
  }
}
