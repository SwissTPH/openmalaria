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
#include "inputData.h"

// submodels:
#include "PkPd/HoshenPkPdModel.h"
#include "PkPd/LSTMPkPdModel.h"

#include <assert.h>

namespace OM { namespace PkPd {
    
PkPdModel::ActiveModel PkPdModel::activeModel = PkPdModel::NON_PKPD;

// weight proportions, used by drug code
const double PkPdModel::wtprop[WithinHost::WithinHostModel::nages] = { 0.116547265, 0.152531009, 0.181214575, 0.202146126, 0.217216287, 0.237405732, 0.257016899, 0.279053187, 0.293361286, 0.309949502, 0.334474135, 0.350044993, 0.371144279, 0.389814144, 0.412366341, 0.453, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };

// -----  static functions  -----

void PkPdModel::init () {
    if (util::ModelOptions::option (util::INCLUDES_PK_PD)) {
	ProteomeManager::init ();
	if (InputData.getScenario().getDrugDescription().present()) {
	    activeModel = LSTM_PKPD;
	    LSTMDrugType::init();
	    LSTMDrug::init ();
	    LSTMPkPdModel::init();
	} else {
	    activeModel = HOSHEN_PKPD;
	    HoshenDrugType::init();
	    HoshenDrug::init ();
	    HoshenPkPdModel::init();
	}
    }
}
void PkPdModel::cleanup () {
    if (activeModel != NON_PKPD) {
	ProteomeManager::cleanup ();
    }
}

PkPdModel* PkPdModel::createPkPdModel () {
    if (activeModel == NON_PKPD) {
	return new PkPdModel();
    } else if (activeModel == LSTM_PKPD) {
	return new LSTMPkPdModel ();
    } else if (activeModel == HOSHEN_PKPD) {
	return new HoshenPkPdModel ();
    }
    assert(false);	// execution shouldn't reach this point
}


double PkPdModel::ageToWeight (double ageYears) {
    return 120.0 * wtprop[WithinHost::WithinHostModel::getAgeGroup(ageYears)];
}

} }