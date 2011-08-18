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
#include "Transmission/PerHost.h"
#include "Transmission/VectorModel.h"
#include "Transmission/Vector/PerHost.h"
#include "inputData.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {
using namespace OM::util;

// -----  PerHost static  -----

AgeGroupInterpolation* PerHost::relAvailAge = AgeGroupInterpolation::dummyObject();
shared_ptr<DecayFunction> PerHost::IRSDecay;
shared_ptr<DecayFunction> PerHost::VADecay;

void PerHost::init () {
    relAvailAge = AgeGroupInterpolation::makeObject( InputData().getModel().getHuman().getAvailabilityToMosquitoes(), "availabilityToMosquitoes" );
}
void PerHost::cleanup () {
    AgeGroupInterpolation::freeObject( relAvailAge );
}

void PerHost::setIRSDescription (const scnXml::IRS& elt) {
    IRSDecay = DecayFunction::makeObject( elt.getDecay(), "IRSDecay" );
}
void PerHost::setVADescription (const scnXml::VectorDeterrent& elt) {
    VADecay = DecayFunction::makeObject( elt.getDecay(), "VADecay" );
}

// -----  PerHost non-static -----

PerHost::PerHost () :
        outsideTransmission(false),
        timestepIRS(TimeStep::never),
        timestepVA(TimeStep::never)
{
    if ( IRSDecay.get() != 0 )
        hetSampleIRS = IRSDecay->hetSample();
    if ( VADecay.get() != 0 )
        hetSampleVA = VADecay->hetSample();
}
void PerHost::initialise (TransmissionModel& tm, double availabilityFactor) {
    _relativeAvailabilityHet = availabilityFactor;
    VectorModel* vTM = dynamic_cast<VectorModel*> (&tm);
    if (vTM != 0) {
        species.resize (vTM->numSpecies);
        for (size_t i = 0; i < vTM->numSpecies; ++i)
            species[i].initialise (vTM->species[i].getHumanBaseParams(), availabilityFactor);
    }
}

void PerHost::setupITN (const TransmissionModel& tm) {
    const VectorModel* vTM = dynamic_cast<const VectorModel*> (&tm);
    if (vTM != 0) {
        net.deploy(vTM->getITNParams());
    }
}
void PerHost::setupIRS () {
    if( IRSDecay.get() == 0 ){
        throw util::xml_scenario_error ("IRS intervention without description of decay");
    }
    timestepIRS = TimeStep::simulation;
}
void PerHost::setupVA () {
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
double PerHost::entoAvailabilityHetVecItv (const Vector::PerHostBase& base, size_t speciesIndex) const {
    double alpha_i = species[speciesIndex].getEntoAvailability();
    if (net.timeOfDeployment() >= TimeStep(0)) {
        alpha_i *= net.relativeAttractiveness(base.net);
    }
    if (timestepIRS >= TimeStep(0))
        alpha_i *= (1.0 - base.IRSDeterrency * IRSDecay->eval (TimeStep::simulation - timestepIRS, hetSampleIRS));
    if (timestepVA >= TimeStep(0))
        alpha_i *= (1.0 - base.VADeterrency * VADecay->eval (TimeStep::simulation - timestepVA, hetSampleVA));

    return alpha_i;
}
double PerHost::probMosqBiting (const Vector::PerHostBase& base, size_t speciesIndex) const {
    double P_B_i = species[speciesIndex].getProbMosqBiting();
    if (net.timeOfDeployment() >= TimeStep(0)) {
        P_B_i *= net.preprandialSurvivalFactor(base.net);
    }
    return P_B_i;
}
double PerHost::probMosqResting (const Vector::PerHostBase& base, size_t speciesIndex) const {
    double pRest = species[speciesIndex].getProbMosqRest();
    if (net.timeOfDeployment() >= TimeStep(0)) {
        pRest *= net.postprandialSurvivalFactor(base.net);
    }
    if (timestepIRS >= TimeStep(0))
        pRest *= (1.0 - base.IRSKillingEffect * IRSDecay->eval (TimeStep::simulation - timestepIRS, hetSampleIRS));
    return pRest;
}

}
}
