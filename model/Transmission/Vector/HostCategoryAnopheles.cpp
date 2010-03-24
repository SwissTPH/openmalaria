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
#include "util/errors.hpp"

namespace OM { namespace Transmission {

void HostCategoryAnopheles::setEntoAvailability(double entoAvailability)
{
	this->entoAvailability = entoAvailability;
}

void HostCategoryAnopheles::setInterventionDescription (const scnXml::Anopheles1& intervDesc) {
  if (intervDesc.getITNDescription().present()) {
    const scnXml::ITNDescription& itnDesc = intervDesc.getITNDescription().get();
    ITNDeterrency = itnDesc.getDeterrency ();
    ITNPreprandialKillingEffect = itnDesc.getPreprandialKillingEffect ();
    ITNPostprandialKillingEffect = itnDesc.getPostprandialKillingEffect ();
  }
  if (intervDesc.getIRSDescription().present()) {
    const scnXml::IRSDescription& irsDesc = intervDesc.getIRSDescription().get();
    IRSDeterrency = irsDesc.getDeterrency ();
    IRSKillingEffect = irsDesc.getKillingEffect ();
  }
  if (intervDesc.getVADescription().present()) {
    const scnXml::VADescription& vaDesc = intervDesc.getVADescription().get();
    VADeterrency = vaDesc.getDeterrency ();
  }
}
void HostCategoryAnopheles::checkInterventionDescriptions (string species) {
    if (InputData.getActiveInterventions()[Interventions::ITN]) {
	if (ITNDeterrency.notSet() ||
	    ITNPreprandialKillingEffect.notSet() ||
	    ITNPostprandialKillingEffect.notSet())
	    throw util::xml_scenario_error (string("ITN intervention without description for ").append(species));
    }
    if (InputData.getActiveInterventions()[Interventions::IRS]) {
	if (IRSDeterrency.notSet() || IRSKillingEffect.notSet())
	    throw util::xml_scenario_error (string("IRS intervention without description for ").append(species));
    }
    if (InputData.getActiveInterventions()[Interventions::VEC_AVAIL]) {
	if (VADeterrency.notSet())
	    throw util::xml_scenario_error (string("Vector Availability intervention without description for ").append(species));
    }
}

} }
