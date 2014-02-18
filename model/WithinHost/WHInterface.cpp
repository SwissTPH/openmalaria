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

#include "WithinHost/WHInterface.h"
#include "WithinHost/WHVivax.h"
#include "WithinHost/DescriptiveWithinHost.h"
#include "WithinHost/DescriptiveIPTWithinHost.h"
#include "WithinHost/CommonWithinHost.h"
#include "WithinHost/Infection/DummyInfection.h"
#include "WithinHost/Infection/EmpiricalInfection.h"
#include "WithinHost/Infection/MolineauxInfection.h"
#include "WithinHost/Infection/PennyInfection.h"
#include "inputData.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
//using namespace std;

#include <cmath>
#include <boost/format.hpp>


namespace OM {
namespace WithinHost {

using namespace OM::util;

// -----  static functions  -----

void WHInterface::init() {
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        WHVivax::init();
    }else{
        WHFalciparum::init();
    }
}

WHInterface* WHInterface::createWithinHostModel () {
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        return new WHVivax();
    }else if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL) ||
            util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL) ||
            util::ModelOptions::option (util::MOLINEAUX_WITHIN_HOST_MODEL) ||
            util::ModelOptions::option (util::PENNY_WITHIN_HOST_MODEL)) {
        return new CommonWithinHost();
    } else if ( util::ModelOptions::option( IPTI_SP_MODEL ) ){
        return new DescriptiveIPTWithinHost();
    }else{
        return new DescriptiveWithinHostModel();
    }
}


// -----  Non-static  -----

WHInterface::WHInterface () :
    numInfs(0)
{
}

WHInterface::~WHInterface()
{
}


void WHInterface::clearInfections (bool) {
    effectiveTreatment();
}

void WHInterface::medicate(string drugAbbrev, double qty, double time, double duration, double bodyMass){
    throw TRACED_EXCEPTION( "should not call medicate() except with CommonWithinHost model", util::Error::WHFeatures );
}

double WHInterface::getTotalDensity() const{
    throw TRACED_EXCEPTION( "should not call getTotalDensity() with non-falciparum model", util::Error::WHFeatures );
}

void WHInterface::continuousIPT (Monitoring::AgeGroup, bool) {
    throw util::xml_scenario_error (string ("Continuous IPT treatment when no IPT description is present in interventions"));
}
void WHInterface::timedIPT (Monitoring::AgeGroup, bool) {
    throw util::xml_scenario_error (string ("Timed IPT treatment when no IPT description is present in interventions"));
}
bool WHInterface::hasIPTiProtection (TimeStep maxInterventionAge) const {
    throw util::xml_scenario_error (string ("Timed IPT treatment when no IPT description is present in interventions"));
}

double WHInterface::getCumulativeh() const{
    throw TRACED_EXCEPTION( "should not call getCumulativeh() with non-falciparum model", util::Error::WHFeatures );
}
double WHInterface::getCumulativeY() const{
    throw TRACED_EXCEPTION( "should not call getCumulativeY() with non-falciparum model", util::Error::WHFeatures );
}


void WHInterface::checkpoint (istream& stream) {
    numInfs & stream;

    if (numInfs > MAX_INFECTIONS)
        throw util::checkpoint_error( (boost::format("numInfs: %1%") %numInfs).str() );
}
void WHInterface::checkpoint (ostream& stream) {
    numInfs & stream;
}

}
}
