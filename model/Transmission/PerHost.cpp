/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2013 Swiss Tropical and Public Health Institute 
 * Copyright (C) 2005-2013 Liverpool School Of Tropical Medicine
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
#include "Transmission/PerHost.h"
#include "Transmission/VectorModel.h"
#include "Transmission/Anopheles/PerHost.h"
#include "util/errors.h"
#include "util/checkpoint.h"

namespace OM {
namespace Transmission {
using namespace OM::util;

// -----  PerHost static  -----

AgeGroupInterpolation* PerHost::relAvailAge = AgeGroupInterpolation::dummyObject();

void PerHost::init ( const scnXml::AgeGroupValues& availabilityToMosquitoes ) {
    relAvailAge = AgeGroupInterpolation::makeObject( availabilityToMosquitoes, "availabilityToMosquitoes" );
}
void PerHost::cleanup () {
    AgeGroupInterpolation::freeObject( relAvailAge );
}

// -----  PerHost non-static -----

PerHost::PerHost (const Transmission::TransmissionModel& tm) :
        outsideTransmission(false),
        net(tm),
        irs(tm)
{
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
void PerHost::setupIRS (const TransmissionModel& tm) {
    const VectorModel* vTM = dynamic_cast<const VectorModel*> (&tm);
    if (vTM != 0) {
        irs.deploy(vTM->getIRSParams());
    }
}

void PerHost::deployEffect( const HumanVectorInterventionEffect& params ){
    //TODO: this adds intervention deployment details; should we ever remove
    // them (once deployment effect is zero or some such)?
    for( ListActiveEffects::iterator it = activeEffects.begin(); it != activeEffects.end(); ++it ){
        if( it->getIndex() == params.getIndex() ){
            // already have a deployment for that description; just update it
            it->redeploy( params );
            return;
        }
    }
    // no deployment for that description: must make a new one
    activeEffects.push_back( params.makeHumanPart() );
}


// Note: in the case an intervention is not present, we can use the approximation
// of Weibull decay over (TimeStep::simulation - TimeStep::never) timesteps
// (easily large enough for conceivable Weibull params that the value is 0.0 when
// rounded to a double. Performance-wise it's perhaps slightly slower than using
// an if() when interventions aren't present.
double PerHost::entoAvailabilityHetVecItv (const Anopheles::PerHostBase& base, size_t speciesIndex) const {
    double alpha_i = species[speciesIndex].getEntoAvailability();
    if (net.timeOfDeployment() >= TimeStep(0)) {
        alpha_i *= net.relativeAttractiveness(base.net);
    }
    if (irs.timeOfDeployment() >= TimeStep(0)) {
        alpha_i *= irs.relativeAttractiveness(base.irs);
    }
    for( ListActiveEffects::const_iterator it = activeEffects.begin(); it != activeEffects.end(); ++it ){
        alpha_i *= it->relativeAttractiveness( speciesIndex );
    }
    return alpha_i;
}
double PerHost::probMosqBiting (const Anopheles::PerHostBase& base, size_t speciesIndex) const {
    double P_B_i = species[speciesIndex].getProbMosqBiting();
    if (net.timeOfDeployment() >= TimeStep(0)) {
        P_B_i *= net.preprandialSurvivalFactor(base.net);
    }
    if (irs.timeOfDeployment() >= TimeStep(0)) {
        P_B_i *= irs.preprandialSurvivalFactor(base.irs);
    }
    for( ListActiveEffects::const_iterator it = activeEffects.begin(); it != activeEffects.end(); ++it ){
        P_B_i *= it->preprandialSurvivalFactor( speciesIndex );
    }
    return P_B_i;
}
double PerHost::probMosqResting (const Anopheles::PerHostBase& base, size_t speciesIndex) const {
    double pRest = species[speciesIndex].getProbMosqRest();
    if (net.timeOfDeployment() >= TimeStep(0)) {
        pRest *= net.postprandialSurvivalFactor(base.net);
    }
    if (irs.timeOfDeployment() >= TimeStep(0)) {
        pRest *= irs.postprandialSurvivalFactor(base.irs);
    }
    for( ListActiveEffects::const_iterator it = activeEffects.begin(); it != activeEffects.end(); ++it ){
        pRest *= it->postprandialSurvivalFactor( speciesIndex );
    }
    return pRest;
}

void PerHost::checkpointIntervs( ostream& stream ){
    activeEffects.size() & stream;
    for( boost::ptr_list<PerHostInterventionData>::iterator it = activeEffects.begin(); it != activeEffects.end(); ++it ){
        *it & stream;
    }
}
void PerHost::checkpointIntervs( istream& stream ){
    size_t l;
    l & stream;
    validateListSize(l);
    activeEffects.clear();
    for( size_t i = 0; i < l; ++i ){
        size_t index;
        index & stream;
        try{
            const interventions::HumanInterventionEffect& gen_params = interventions::InterventionManager::instance->getEffect( index );   // may throw
            const HumanVectorInterventionEffect *params = dynamic_cast<const HumanVectorInterventionEffect*>( &gen_params );
            if( params == 0 )
                throw util::base_exception( "" );       // see catch block below
            PerHostInterventionData *v = params->makeHumanPart( stream, index );
            activeEffects.push_back( v );
        }catch( util::base_exception e ){
            // two causes, both boil down to index being wrong
            throw util::checkpoint_error( "bad value in checkpoint file" );
        }
    }
}

}
}
