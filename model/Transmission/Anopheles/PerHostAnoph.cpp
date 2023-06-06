/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
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

#include "Transmission/Anopheles/PerHostAnoph.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {
namespace Anopheles {

vector<PerHostAnophParams> PerHostAnophParams::params;

// ----- Per host, per species, non-static -----

PerHostAnophParams::PerHostAnophParams (const scnXml::Mosq& mosq) {
    const string &distr = mosq.getAvailability().getDistr();
    if(distr == "const" || distr == "lognormal")
        entoAvailability = make_unique<util::LognormalSampler>(1.0, mosq.getAvailability());
    else if(distr == "gamma")
        entoAvailability = make_unique<util::GammaSampler>(mosq.getAvailability());
    else
        throw util::xml_scenario_error( "error ento availability: unknown distirbution "+distr);

    probMosqBiting.setParams( mosq.getMosqProbBiting() );
    probMosqFindRestSite.setParams( mosq.getMosqProbFindRestSite() );
    probMosqSurvivalResting.setParams( mosq.getMosqProbResting() );
}

void PerHostAnoph::initialise (LocalRng& rng, size_t species, double availabilityFactor)
{
    const PerHostAnophParams& base = PerHostAnophParams::get(species);
    entoAvailability = base.entoAvailability->sample(rng) * availabilityFactor;
    probMosqBiting = base.probMosqBiting.sample(rng);
    auto pRest1 = base.probMosqFindRestSite.sample(rng);
    auto pRest2 = base.probMosqSurvivalResting.sample(rng);
    probMosqRest = pRest1 * pRest2;
}

}
}
}
