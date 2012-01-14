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
#include "Transmission/PerHostTransmission.h"
#include "Transmission/Vector/VectorTransmission.h"
#include "Transmission/Vector/AnophelesHumanParams.h"
#include "inputData.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {
using namespace OM::util;

// -----  PerHostTransmission static  -----

AgeGroupInterpolation* PerHostTransmission::relAvailAge = AgeGroupInterpolation::dummyObject();
shared_ptr<DecayFunction> PerHostTransmission::VADecay;

void PerHostTransmission::init () {
    relAvailAge = AgeGroupInterpolation::makeObject( InputData().getModel().getHuman().getAvailabilityToMosquitoes(), "availabilityToMosquitoes" );
}
void PerHostTransmission::cleanup () {
    AgeGroupInterpolation::freeObject( relAvailAge );
}

void PerHostTransmission::setVADescription (const scnXml::VectorDeterrent& elt) {
    VADecay = DecayFunction::makeObject( elt.getDecay(), "VADecay" );
}

// -----  PerHostTransmission non-static -----

PerHostTransmission::PerHostTransmission (const Transmission::TransmissionModel& tm) :
        outsideTransmission(false),
        timestepVA(TimeStep::never),
        net(tm),
        irs(tm)
{
    if ( VADecay.get() != 0 )
        hetSampleVA = VADecay->hetSample();
}
void PerHostTransmission::initialise (TransmissionModel& tm, double availabilityFactor) {
    _relativeAvailabilityHet = availabilityFactor;
    VectorTransmission* vTM = dynamic_cast<VectorTransmission*> (&tm);
    if (vTM != 0) {
        species.resize (vTM->numSpecies);
        for (size_t i = 0; i < vTM->numSpecies; ++i)
            species[i].initialise (vTM->species[i].getHumanBaseParams(), availabilityFactor);
    }
}

void PerHostTransmission::setupITN (const TransmissionModel& tm) {
    const VectorTransmission* vTM = dynamic_cast<const VectorTransmission*> (&tm);
    if (vTM != 0) {
        net.deploy(vTM->getITNParams());
    }
}
void PerHostTransmission::setupIRS (const TransmissionModel& tm) {
    const VectorTransmission* vTM = dynamic_cast<const VectorTransmission*> (&tm);
    if (vTM != 0) {
        irs.deploy(vTM->getIRSParams());
    }
}
void PerHostTransmission::setupVA () {
    if( VADecay.get() == 0 ){
        throw util::xml_scenario_error ("Vector availability intervention without description of decay");
    }
    timestepVA = TimeStep::simulation;
}


// Note: in the case an intervention is not present, we can use the approximation
// of Weibull decay over (TimeStep::simulation - TimeStep::never) timesteps
// (easily large enough for conceivable Weibull params that the value is 0.0 when
// rounded to a double. Performance-wise it's perhaps slightly slower than using
// an if() when interventions aren't present.
double PerHostTransmission::entoAvailabilityHetVecItv (const AnophelesHumanParams& base, size_t speciesIndex) const {
    double alpha_i = species[speciesIndex].entoAvailability;
    if (net.timeOfDeployment() >= TimeStep(0)) {
        alpha_i *= net.relativeAttractiveness(base.net);
    }
    if (irs.timeOfDeployment() >= TimeStep(0)) {
        alpha_i *= irs.relativeAttractiveness(base.irs);
    }
    if (timestepVA >= TimeStep(0)) {
        alpha_i *= (1.0 - base.VADeterrency * VADecay->eval (TimeStep::simulation - timestepVA, hetSampleVA));
    }
    return alpha_i;
}
double PerHostTransmission::probMosqBiting (const AnophelesHumanParams& base, size_t speciesIndex) const {
    double P_B_i = species[speciesIndex].probMosqBiting;
    if (net.timeOfDeployment() >= TimeStep(0)) {
        P_B_i *= net.preprandialSurvivalFactor(base.net);
    }
    if (irs.timeOfDeployment() >= TimeStep(0)) {
        P_B_i *= irs.preprandialSurvivalFactor(base.irs);
    }
    return P_B_i;
}
double PerHostTransmission::probMosqResting (const AnophelesHumanParams& base, size_t speciesIndex) const {
    double pRest = species[speciesIndex].probMosqRest;
    if (net.timeOfDeployment() >= TimeStep(0)) {
        pRest *= net.postprandialSurvivalFactor(base.net);
    }
    if (irs.timeOfDeployment() >= TimeStep(0)) {
        pRest *= irs.postprandialSurvivalFactor(base.irs);
    }
    return pRest;
}


// ----- HostMosquitoInteraction non-static -----

void HostMosquitoInteraction::initialise (const AnophelesHumanParams& base, double availabilityFactor)
{
    entoAvailability = base.entoAvailability.sample() * availabilityFactor;
    probMosqBiting = base.probMosqBiting.sample();
    probMosqRest = base.probMosqFindRestSite.sample() * base.probMosqSurvivalResting.sample();
}

}
}
