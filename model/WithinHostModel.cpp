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

#include "WithinHostModel.h"
#include "WithinHostModel/Descriptive.h"
#include "WithinHostModel/OldIPT.h"
#include "WithinHostModel/Dummy.h"
#include "WithinHostModel/Empirical.h"
#include "inputData.h"
#include "WithinHostModel/DescriptiveInfection.h"
#include <stdexcept>

using namespace std;

// weight proportions, used by drug code
const double WithinHostModel::wtprop[nwtgrps] = { 0.116547265, 0.152531009, 0.181214575, 0.202146126, 0.217216287, 0.237405732, 0.257016899, 0.279053187, 0.293361286, 0.309949502, 0.334474135, 0.350044993, 0.371144279, 0.389814144, 0.412366341, 0.453, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
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
    EmpiricalInfection::initParameters();
  } else {
    DescriptiveInfection::initParameters ();
    OldIPTWithinHostModel::initParameters();
  }
}

void WithinHostModel::clear() {
  OldIPTWithinHostModel::clearParameters();
  DescriptiveInfection::clearParameters();
}

WithinHostModel* WithinHostModel::createWithinHostModel () {
  if (Global::modelVersion & DUMMY_WITHIN_HOST_MODEL) {
    return new DummyWithinHostModel();
  } else if (Global::modelVersion & EMPIRICAL_WITHIN_HOST_MODEL) {
    return new EmpiricalWithinHostModel();
  } else {
    if (OldIPTWithinHostModel::iptActive)
      return new OldIPTWithinHostModel();
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
    if (OldIPTWithinHostModel::iptActive)
      return new OldIPTWithinHostModel(in);
    else
      return new DescriptiveWithinHostModel(in);
  }
}

WithinHostModel::WithinHostModel(istream& in) {
  in >> _cumulativeInfections; 
  in >> _pTransToMosq; 
  in >> totalDensity;
  in >> timeStepMaxDensity;
}

void WithinHostModel::clearInfections (bool) {
  clearAllInfections();
}

void WithinHostModel::IPTiTreatment (double compliance, int ageGroup) {
  throw xml_scenario_error (string ("Timed IPT treatment when no IPT description is present in interventions"));
}
