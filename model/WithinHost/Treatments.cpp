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

#include "WithinHost/Treatments.h"
#include "util/errors.h"
#include "util/random.h"
#include "schema/healthSystem.h"

namespace OM {
namespace WithinHost {

// ———  static  ———

TreatmentId TreatmentId::legacy;

boost::ptr_vector<vector<Treatments> > Treatments::groups;

void Treatments::init(){
    if( TimeStep::interval == 5 ){
        vector<Treatments>* legacyGroup = new vector<Treatments>();
        Treatments treat( 1.0 );
        treat.effects.push_back( Action( -1, BOTH ) );
        legacyGroup->push_back( treat );
        TreatmentId::legacy.id = groups.size();  // i.e. the index used below; probably 0
        groups.push_back( legacyGroup );
    }else{
        TreatmentId::legacy.id = numeric_limits<uint32_t>::max();
    }
}

TreatmentId Treatments::addTreatment(const scnXml::TreatmentDescription& desc){
    vector<Treatments>* group = new vector<Treatments>();
    group->reserve( desc.getOption().size() );
    
    double cumP = 0.0;
    for( scnXml::TreatmentDescription::OptionConstIterator it =
        desc.getOption().begin(), end = desc.getOption().end(); it != end; ++it )
    {
        cumP += it->getPSelection();
        Treatments treat( cumP, *it );
        group->push_back( treat );
    }
    
    // we expect the prob. to be roughly one as an error check, but allow slight deviation
    if( cumP < 0.99 || cumP > 1.01 ) throw util::xml_scenario_error( "sum of pSelection of a group of treatments is not 1" );
    for( vector<Treatments>::iterator it = group->begin(), end = group->end(); it != end; ++it ){
        it->pCum /= cumP;
    }
    
    TreatmentId id;
    id.id = groups.size();
    groups.push_back( group );
    return id;
}

const Treatments& Treatments::select(TreatmentId treatId){
    assert( treatId.id >= 0 && treatId.id < groups.size() );
    vector<Treatments>& group = groups[treatId.id];
    assert( group.size() >= 1 );
    if( group.size() == 1 ) return group[0];        // take only option
    
    double x = util::random::uniform_01();      // random sample: choose
    for( vector<Treatments>::const_iterator it = group.begin(), end =
        group.end(); it != end; ++it ){
        if( it->pCum > x ) return *it;
    }
    assert( false );    // last item should have pCum=1 and x<1 in theory
    return group[0]; // remain type safe
}


// ———   non-static  ———

Treatments::Treatments( double pCum ) :
    pCum( pCum )
{}

Treatments::Stages stageFromString( const std::string& str ){
    if( str == "liver" ) return Treatments::LIVER;
    if( str == "blood" ) return Treatments::BLOOD;
    assert( str == "both" );
    return Treatments::BOTH;
}
Treatments::Treatments( double pCum, const scnXml::TreatmentOption& elt ) :
    pCum( pCum )
{
    for( scnXml::TreatmentOption::TreatmentConstIterator it =
        elt.getTreatment().begin(), end = elt.getTreatment().end(); it != end;
        ++it )
    {
        int len = it->getTimesteps();
        if( len < -1 || len == 0 ){ 
            throw util::xml_scenario_error( "prophylacticTreatment: timesteps: must be ≥ 1 or have special value -1" );
        }
        Stages stage = stageFromString( it->getStage() );
        if( TimeStep::interval != 5 && stage != BOTH ){
            throw util::unimplemented_exception(
                "differentiation of infection stages for simple treatment (alternative: use the PK/PD model)" );
        }
        effects.push_back( Action( len, stage ) );
    }
}

}
}
