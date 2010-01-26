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

#include "PkPd/PkPdModel.h"
#include "PkPd/Drug/HoshenDrugType.h"
#include "Global.h"
#include "PkPd/Proteome.h"
#include "util/ModelOptions.hpp"

// submodels:
#include "PkPd/HoshenPkPdModel.h"
#include "PkPd/LSTMPkPdModel.h"

namespace OM { namespace PkPd {
    
// Temporary switch to use the LSTM model − this may eventually be determined by XML data or XML model version.
const bool Use_LSTM = false;

// weight proportions, used by drug code
const double PkPdModel::wtprop[WithinHost::WithinHostModel::nages] = { 0.116547265, 0.152531009, 0.181214575, 0.202146126, 0.217216287, 0.237405732, 0.257016899, 0.279053187, 0.293361286, 0.309949502, 0.334474135, 0.350044993, 0.371144279, 0.389814144, 0.412366341, 0.453, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };

// -----  static functions  -----

void PkPdModel::init () {
  if (util::ModelOptions::option (util::INCLUDES_PK_PD)) {
    ProteomeManager::init ();
    if (Use_LSTM) {
      LSTMDrugType::init();
      LSTMDrug::init ();
      LSTMPkPdModel::init();
    }
    else {
      HoshenDrugType::init();
      HoshenDrug::init ();
      HoshenPkPdModel::init();
    }
  }
}
void PkPdModel::cleanup () {
    if (util::ModelOptions::option (util::INCLUDES_PK_PD)) {
	ProteomeManager::cleanup ();
    }
}

void PkPdModel::staticCheckpoint (istream& stream) {
  if (util::ModelOptions::option (util::INCLUDES_PK_PD)) {
//     ProteomeManager::read (in);
  }
}
void PkPdModel::staticCheckpoint (ostream& stream) {
  if (util::ModelOptions::option (util::INCLUDES_PK_PD)) {
//     ProteomeManager::write (out);
  }
}

PkPdModel* PkPdModel::createPkPdModel () {
  if (util::ModelOptions::option (util::INCLUDES_PK_PD)) {
    if (Use_LSTM)
      return new LSTMPkPdModel ();
    else
      return new HoshenPkPdModel ();
  }
  return new PkPdModel();
}


double PkPdModel::ageToWeight (double ageYears) {
    return 120.0 * wtprop[WithinHost::WithinHostModel::getAgeGroup(ageYears)];
}

} }