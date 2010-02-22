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

#include "WithinHost/DummyInfection.h"
#include "WithinHost/Common.h"
#include "inputData.h"
#include "util/gsl.h"
#include "util/ModelOptions.hpp"

#include <algorithm>
#include <sstream>
#include <string.h>

namespace OM { namespace WithinHost {
    
CommonInfection* createDummyInfection (uint32_t protID) {
    return new DummyInfection (protID);
}
CommonInfection* checkpointedDummyInfection (istream& stream) {
    return new DummyInfection (stream);
}

void DummyInfection::initParameters () {
    CommonWithinHost::createInfection = &createDummyInfection;
    CommonWithinHost::checkpointedInfection = &checkpointedDummyInfection;
}

DummyInfection::DummyInfection (uint32_t protID) :
    CommonInfection(protID)
{
    _density=16;	// increased by DH to avoid zeros in initialKappa
}

bool DummyInfection::updateDensity(int simulationTime, double survivalFactor) {
    const double GROWTH_RATE = 8.0;
    const double PARASITE_THRESHOLD = 1;
    
    _density = (int(_density*GROWTH_RATE) % 20000) * survivalFactor;
    _cumulativeExposureJ += Global::interval * _density;
    
    if (_density < PARASITE_THRESHOLD) {
	return true;
    } else {
	return false;
    }
}


DummyInfection::DummyInfection (istream& stream) :
    CommonInfection (stream)
{}

} }