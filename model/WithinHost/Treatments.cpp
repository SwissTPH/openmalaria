/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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
#include "util/timeConversions.h"
#include "schema/healthSystem.h"

namespace OM {
namespace WithinHost {

// ———  static  ———

boost::ptr_vector<Treatments> Treatments::treatments;

TreatmentId Treatments::addTreatment( const scnXml::TreatmentOption& desc ){
    Treatments* treatment = new Treatments( desc );
    TreatmentId id( treatments.size() );
    treatments.push_back( treatment );
    return id;
}


// ———   non-static  ———

Treatments::Stages stageFromString( const std::string& str ){
    if( str == "liver" ) return Treatments::LIVER;
    if( str == "blood" ) return Treatments::BLOOD;
    if( str == "both" ) return Treatments::BOTH;
    throw util::format_error( std::string("stage must be liver, blood or both, not ").append(str) );
}
Treatments::Treatments( const scnXml::TreatmentOption& elt ) :
    TriggeredDeployments(elt), timeLiver(SimTime::zero()), timeBlood(SimTime::zero())
{
    for( scnXml::TreatmentOption::ClearInfectionsConstIterator it =
        elt.getClearInfections().begin(), end = elt.getClearInfections().end();
        it != end; ++it )
    {
        try{
            //NOTE: if changing XSD, this should not be called "timesteps" or have a default unit
            SimTime len = UnitParse::readShortDuration( it->getTimesteps(), UnitParse::STEPS );
            if( len < -SimTime::oneTS() || len == SimTime::zero() ){ 
                throw util::format_error( "timesteps must be ≥ 1 or have special value -1" );
            }
            Stages stage = stageFromString( it->getStage() );
            if( stage & LIVER ){
                if( timeLiver != SimTime::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of liver stage effect" );
                timeLiver = len;
            }
            if( stage & BLOOD ){
                if( timeBlood != SimTime::zero() )   // existing treatment configuration
                    throw util::format_error( "multiple specification of blood stage effect" );
                timeBlood = len;
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string(".../clearInfections/timesteps: ").append(e.message()) );
        }
    }
}

}
}
