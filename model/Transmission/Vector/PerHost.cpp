/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Transmission/Vector/PerHost.h"
#include "inputData.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {
namespace Vector {

// ----- PerHostPerSpecies non-static -----

void PerHost::initialise (const PerHostBase& base, double availabilityFactor)
{
    entoAvailability = base.entoAvailability.sample() * availabilityFactor;
    probMosqBiting = base.probMosqBiting.sample();
    probMosqRest = base.probMosqFindRestSite.sample() * base.probMosqSurvivalResting.sample();
}

void PerHostBase::operator =(const scnXml::Mosq& mosq)
{
    entoAvailability.setParams( numeric_limits<double>::quiet_NaN(),
                                mosq.getAvailabilityVariance().getValue() );
    probMosqBiting.setParams( mosq.getMosqProbBiting() );
    probMosqFindRestSite.setParams( mosq.getMosqProbFindRestSite() );
    probMosqSurvivalResting.setParams( mosq.getMosqProbResting() );
}

void PerHostBase::setITNDescription (const ITNParams& params,
                                     const scnXml::ITNDescription::AnophelesParamsType& elt,
                                     double proportionUse){
    net.init( params, elt, proportionUse );
}

void PerHostBase::setIRSDescription (const scnXml::IRSDescription& irsDesc){
    IRSDeterrency = irsDesc.getDeterrency().getValue();
    IRSKillingEffect = irsDesc.getKillingEffect().getValue();
}

void PerHostBase::setVADescription (const scnXml::BaseInterventionDescription& vaDesc) {
    VADeterrency = vaDesc.getDeterrency().getValue();
}

}
}
}
