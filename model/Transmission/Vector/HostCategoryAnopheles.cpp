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

namespace OM { namespace Transmission {

void HostCategoryAnopheles::setEntoAvailability(double entoAvailability)
{
	this->entoAvailability = entoAvailability;
}

void HostCategoryAnopheles::setInterventionDescription (const scnXml::Anopheles1& intervDesc) {
  if (intervDesc.getITNDescription().present()) {
    const scnXml::ITNDescription& itnDesc = intervDesc.getITNDescription().get();
    ITNDeterrency = util::DecayFunction::makeObject( itnDesc.getDeterrency(), "deterrency" );
    ITNPreprandialKillingEffect = util::DecayFunction::makeObject( itnDesc.getPreprandialKillingEffect(), "preprandialKillingEffect" );
    ITNPostprandialKillingEffect = util::DecayFunction::makeObject( itnDesc.getPostprandialKillingEffect(), "postprandialKillingEffect" );
  }
  if (intervDesc.getIRSDescription().present()) {
    const scnXml::IRSDescription& irsDesc = intervDesc.getIRSDescription().get();
    IRSDeterrency = util::DecayFunction::makeObject( irsDesc.getDeterrency(), "deterrency" );
    IRSKillingEffect = util::DecayFunction::makeObject( irsDesc.getKillingEffect(), "killingEffect" );
  }
  if (intervDesc.getVADescription().present()) {
    const scnXml::BaseInterventionDescription& vaDesc = intervDesc.getVADescription().get();
    VADeterrency = util::DecayFunction::makeObject( vaDesc.getDeterrency(), "deterrency" );
  }
}
void HostCategoryAnopheles::checkInterventionDescriptions (string species) {
    if (InputData.getActiveInterventions()[Interventions::ITN]) {
	if (ITNDeterrency.get()==0 ||
	    ITNPreprandialKillingEffect.get()==0 ||
	    ITNPostprandialKillingEffect.get()==0)
	    throw util::xml_scenario_error (string("ITN intervention without description for ").append(species));
    }
    if (InputData.getActiveInterventions()[Interventions::IRS]) {
	if (IRSDeterrency.get()==0 || IRSKillingEffect.get()==0)
	    throw util::xml_scenario_error (string("IRS intervention without description for ").append(species));
    }
    if (InputData.getActiveInterventions()[Interventions::VEC_AVAIL]) {
	if (VADeterrency.get()==0)
	    throw util::xml_scenario_error (string("Vector Availability intervention without description for ").append(species));
    }
}

} }
