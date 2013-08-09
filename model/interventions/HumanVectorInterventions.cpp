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

#include "interventions/HumanVectorInterventions.h"
 
namespace OM { namespace interventions {

void HumanVectorInterventions::deploy( const HumanVectorInterventionParams& params ){
    //TODO: this adds intervention deployment details; should we ever remove
    // them (once deployment effect is zero or some such)?
    for( ActiveList::iterator it = active.begin(); it != active.end(); ++it ){
        if( it->getIndex() == params.getIndex() ){
            // already have a deployment for that description; just update it
            it->deploy( params );
            return;
        }
    }
    // no deployment for that description: must make a new one
    HumanVectorIntervention *deploy = params.makeHumanPart();
    active.push_back( deploy );
    deploy->deploy( params );
}

double HumanVectorInterventions::relativeAttractiveness( size_t speciesIndex )const{
    double factor = 1;
    for( ActiveList::const_iterator it = active.begin(); it != active.end(); ++it ){
        factor *= it->relativeAttractiveness( InterventionManager::instance->getEffect(it->getIndex()), speciesIndex );
    }
    return factor;
}

double HumanVectorInterventions::preprandialSurvivalFactor( size_t speciesIndex )const{
    double factor = 1;
    for( ActiveList::const_iterator it = active.begin(); it != active.end(); ++it ){
        factor *= it->preprandialSurvivalFactor( InterventionManager::instance->getEffect(it->getIndex()), speciesIndex );
    }
    return factor;
}

double HumanVectorInterventions::postprandialSurvivalFactor( size_t speciesIndex )const{
    double factor = 1;
    for( ActiveList::const_iterator it = active.begin(); it != active.end(); ++it ){
        factor *= it->postprandialSurvivalFactor( InterventionManager::instance->getEffect(it->getIndex()), speciesIndex );
    }
    return factor;
}

void HumanVectorInterventions::checkpoint( ostream& stream ){
    active.size() & stream;
    for( boost::ptr_list<HumanVectorIntervention>::iterator it = active.begin(); it != active.end(); ++it ){
        *it & stream;
    }
}
void HumanVectorInterventions::checkpoint( istream& stream ){
    size_t l;
    l & stream;
    validateListSize(l);
    active.clear();
    for( size_t i = 0; i < l; ++i ){
        size_t index;
        index & stream;
        try{
            const HumanInterventionEffect& gen_params = InterventionManager::instance->getEffect( index );   // may throw
            const HumanVectorInterventionParams *params = dynamic_cast<const HumanVectorInterventionParams*>( &gen_params );
            if( params == 0 )
                throw util::base_exception( "" );       // see catch block below
            HumanVectorIntervention *v = params->makeHumanPart( stream, index );
            active.push_back( v );
        }catch( util::base_exception e ){
            // two causes, both boil down to index being wrong
            throw util::checkpoint_error( "bad value in checkpoint file" );
        }
    }
}

} }
