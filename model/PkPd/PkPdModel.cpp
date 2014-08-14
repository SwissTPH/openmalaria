/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#include "PkPd/PkPdModel.h"
// #include "PkPd/Drug/HoshenDrugType.h"
#include "Global.h"
// #include "PkPd/Proteome.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include <schema/scenario.h>

// submodels:
// #include "PkPd/HoshenPkPdModel.h"
#include "PkPd/LSTMPkPdModel.h"
#include "PkPd/VoidPkPdModel.h"
#include "PkPd/LSTMTreatments.h"

#include <assert.h>
#include <stdexcept>
#include <limits.h>

namespace OM { namespace PkPd {

PkPdModel::ActiveModel PkPdModel::activeModel = PkPdModel::NON_PKPD;


// -----  static functions  -----

void PkPdModel::init( const scnXml::Scenario& scenario ){
    if (util::ModelOptions::option (util::INCLUDES_PK_PD)) {
        if (scenario.getPharmacology().present()) {
            activeModel = LSTM_PKPD;
            LSTMDrugType::init(scenario.getPharmacology().get().getDrugs());
            LSTMTreatments::init(scenario.getPharmacology().get().getTreatments());
        } else {
            throw util::xml_scenario_error( "pharmacology element required in XML" );
        }
        /* if ... {
            // Hoshen model has been removed.
            activeModel = HOSHEN_PKPD;
            ProteomeManager::init ();
            HoshenDrugType::init();
        } */
    }
}

/*
void PkPdModel::cleanup () {
    if (activeModel == HOSHEN_PKPD) {
        assert( false );
        HoshenDrugType::cleanup();
        ProteomeManager::cleanup ();
    }
}
*/

PkPdModel* PkPdModel::createPkPdModel () {
    if (activeModel == NON_PKPD) {
        return new VoidPkPdModel();
    } else if (activeModel == LSTM_PKPD) {
        return new LSTMPkPdModel ();
    } /* else if (activeModel == HOSHEN_PKPD) {
        return new HoshenPkPdModel ();
    } */
    throw TRACED_EXCEPTION_DEFAULT("bad PKPD model");
}


} }
