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
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/random.h"
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
    throw util::xml_scenario_error( std::string("treatment action: stage must be liver, blood or both, not ").append(str) );
}
Treatments::Treatments( const scnXml::TreatmentOption& elt ) :
    TriggeredDeployments(elt), timeLiver(sim::zero()), timeBlood(sim::zero())
{
    for( scnXml::TreatmentOption::ClearInfectionsConstIterator it =
        elt.getClearInfections().begin(), end = elt.getClearInfections().end();
        it != end; ++it )
    {
        //FIXME(schema): input should not be in time steps
        int len = it->getTimesteps();
        if( len < -1 || len == 0 ){ 
            throw util::xml_scenario_error( "prophylacticTreatment: timesteps: must be ≥ 1 or have special value -1" );
        }
        Stages stage = stageFromString( it->getStage() );
        if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
            if( stage != BLOOD || len != -1 )
                throw util::unimplemented_exception( "vivax model requires treatments configured as blood-stage with timesteps=-1" );
            // Actually, the model ignores these parameters; we just don't want somebody thinking it doesn't.
        }
        if( stage & LIVER ){
            if( timeLiver != sim::zero() )   // existing treatment configuration
                throw util::xml_scenario_error( "treatment action: multiple specification of liver stage effect" );
            timeLiver = sim::fromTS(len);
        }
        if( stage & BLOOD ){
            if( timeBlood != sim::zero() )   // existing treatment configuration
                throw util::xml_scenario_error( "treatment action: multiple specification of blood stage effect" );
            timeBlood = sim::fromTS(len);
        }
    }
}

}
}
