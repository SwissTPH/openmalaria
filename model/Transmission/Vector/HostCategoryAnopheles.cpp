/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Transmission/Vector/HostCategoryAnopheles.h"
#include "inputData.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {

void HostCategoryAnopheles::setEntoAvailability(double entoAvailability)
{
    this->entoAvailability = entoAvailability;
}

void HostCategoryAnopheles::setInterventionDescription (const scnXml::Anopheles& intervDesc, const string& species) {
    assert(false);
    /*FIXME
    if ( InputData.isInterventionActive(Interventions::ITN) ) {
        if (!intervDesc.getITNDescription().present()) {
            throw util::xml_scenario_error (string("ITN intervention without description for ").append(species));
        }
        const scnXml::ITNDescription& itnDesc = intervDesc.getITNDescription().get();
        ITNDeterrency = itnDesc.getDeterrency().getValue();
        ITNPreprandialKillingEffect = itnDesc.getPreprandialKillingEffect().getValue();
        ITNPostprandialKillingEffect = itnDesc.getPostprandialKillingEffect().getValue();
    }
    if (InputData.isInterventionActive(Interventions::IRS)) {
        if (!intervDesc.getIRSDescription().present()) {
            throw util::xml_scenario_error (string("IRS intervention without description for ").append(species));
        }
        const scnXml::IRSDescription& irsDesc = intervDesc.getIRSDescription().get();
        IRSDeterrency = irsDesc.getDeterrency().getValue();
        IRSKillingEffect = irsDesc.getKillingEffect().getValue();
    }
    if (InputData.isInterventionActive(Interventions::VEC_AVAIL)) {
        if (!intervDesc.getVADescription().present()) {
            throw util::xml_scenario_error (string("Vector Availability intervention without description for ").append(species));
        }
        const scnXml::BaseInterventionDescription& vaDesc = intervDesc.getVADescription().get();
        VADeterrency = vaDesc.getDeterrency().getValue();
    }
    */
}

}
}
