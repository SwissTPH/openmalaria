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
#include "Transmission/PerHost.h"
#include "Transmission/VectorModel.h"
#include "Transmission/Anopheles/PerHost.h"
#include "interventions/InterventionManager.hpp"
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
        _relativeAvailabilityHet(numeric_limits<double>::signaling_NaN())
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

void PerHost::update(){
    for( ListActiveComponents::iterator it = activeComponents.begin(); it != activeComponents.end(); ++it ){
        it->update();
    }
}

void PerHost::deployComponent( const HumanVectorInterventionComponent& params ){
    // This adds per-host per-intervention details to the host's data set.
    // This data is never removed since it can contain per-host heterogeneity samples.
    for( ListActiveComponents::iterator it = activeComponents.begin(); it != activeComponents.end(); ++it ){
        if( it->id() == params.id() ){
            // already have a deployment for that description; just update it
            it->redeploy( params );
            return;
        }
    }
    // no deployment for that description: must make a new one
    activeComponents.push_back( params.makeHumanPart() );
}


// Note: in the case an intervention is not present, we can use the approximation
// of Weibull decay over (TimeStep::simulation - TimeStep::never) timesteps
// (easily large enough for conceivable Weibull params that the value is 0.0 when
// rounded to a double. Performance-wise it's perhaps slightly slower than using
// an if() when interventions aren't present.
double PerHost::entoAvailabilityHetVecItv (const Anopheles::PerHostBase& base,
                                size_t speciesIndex) const {
    double alpha_i = species[speciesIndex].getEntoAvailability();
    for( ListActiveComponents::const_iterator it = activeComponents.begin(); it != activeComponents.end(); ++it ){
        alpha_i *= it->relativeAttractiveness( speciesIndex );
    }
    return alpha_i;
}
double PerHost::probMosqBiting (const Anopheles::PerHostBase& base, size_t speciesIndex) const {
    double P_B_i = species[speciesIndex].getProbMosqBiting();
    for( ListActiveComponents::const_iterator it = activeComponents.begin(); it != activeComponents.end(); ++it ){
        P_B_i *= it->preprandialSurvivalFactor( speciesIndex );
    }
    return P_B_i;
}
double PerHost::probMosqResting (const Anopheles::PerHostBase& base, size_t speciesIndex) const {
    double pRest = species[speciesIndex].getProbMosqRest();
    for( ListActiveComponents::const_iterator it = activeComponents.begin(); it != activeComponents.end(); ++it ){
        pRest *= it->postprandialSurvivalFactor( speciesIndex );
    }
    return pRest;
}

void PerHost::checkpointIntervs( ostream& stream ){
    activeComponents.size() & stream;
    for( boost::ptr_list<PerHostInterventionData>::iterator it = activeComponents.begin(); it != activeComponents.end(); ++it ){
        *it & stream;
    }
}
void PerHost::checkpointIntervs( istream& stream ){
    size_t l;
    l & stream;
    validateListSize(l);
    activeComponents.clear();
    for( size_t i = 0; i < l; ++i ){
        interventions::ComponentId id( stream );
        try{
            const interventions::HumanInterventionComponent& gen_params =
                interventions::InterventionManager::getComponent( id );   // may throw
            const HumanVectorInterventionComponent *params =
                dynamic_cast<const HumanVectorInterventionComponent*>( &gen_params );
            if( params == 0 )
                throw util::base_exception( "" );       // see catch block below
            PerHostInterventionData *v = params->makeHumanPart( stream, id );
            activeComponents.push_back( v );
        }catch( util::base_exception e ){
            // two causes, both boil down to index being wrong
            throw util::checkpoint_error( "bad value in checkpoint file" );
        }
    }
}

}
}
