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

#include "withinHostModel.h"
#include "oldWithinHostModel.h"
#include "OldIPTWithinHostModel.h"
#include "DummyWithinHostModel.h"
#include "inputData.h"

using namespace std;

// weight proportions, used by drug code
const double WithinHostModel::wtprop[nwtgrps] = { 0.116547265, 0.152531009, 0.181214575, 0.202146126, 0.217216287, 0.237405732, 0.257016899, 0.279053187, 0.293361286, 0.309949502, 0.334474135, 0.350044993, 0.371144279, 0.389814144, 0.412366341, 0.453, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
double WithinHostModel::sigma_i;
double WithinHostModel::immPenalty_22;
double WithinHostModel::asexImmRemain;
double WithinHostModel::immEffectorRemain;

// -----  Initialization  -----

void WithinHostModel::init() {
  sigma_i=sqrt(getParameter(Params::SIGMA_I_SQ));
  immPenalty_22=1-exp(getParameter(Params::IMMUNITY_PENALTY));
  immEffectorRemain=exp(-getParameter(Params::IMMUNE_EFFECTOR_DECAY));
  asexImmRemain=exp(-getParameter(Params::ASEXUAL_IMMUNITY_DECAY));
}

WithinHostModel* WithinHostModel::createWithinHostModel () {
  if (Global::modelVersion & DUMMY_WITHIN_HOST_MODEL) {
    return new DummyWithinHostModel();
  } else {
    if (IPTIntervention::IPT)
      return new OldIPTWithinHostModel();
    else
      return new OldWithinHostModel();
  }
}

void WithinHostModel::IPTClearInfections (Event&) {
  clearAllInfections();
}
