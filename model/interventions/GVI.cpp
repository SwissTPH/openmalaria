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

#include "interventions//GVI.h"
//TODO: we shouldn't have a dependency on the vector/transmission model class
//here; currently it's a work-around for GVI parameters not always being present.
#include "Transmission/VectorModel.h"
#include "util/random.h"
#include "util/errors.h"
#include "R_nmath/qnorm.h"
#include <cmath>

namespace OM { namespace interventions {
    using util::random::poisson;

void GVIParams::init( const scnXml::GVIDescription& elt) {
    decay = DecayFunction::makeObject( elt.getDecay(), "interventions.human.vector.decay" );
}

void GVIAnophelesParams::init(
    const GVIParams& params,
    const scnXml::GVIDescription::AnophelesParamsType& elt)
{
    if( _relativeAttractiveness == _relativeAttractiveness ){
        throw util::unimplemented_exception( "multiple GVI interventions" );
    }
    _relativeAttractiveness = elt.getDeterrency().getValue();
    _preprandialKillingEffect = elt.getPreprandialKillingEffect().getValue();
    _postprandialKillingEffect = elt.getPostprandialKillingEffect().getValue();
    // Simpler version of ITN usage/action:
    double propActive = elt.getPropActive();
    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = propActive;
    proportionUnprotected = 1.0 - proportionProtected;
}


// ———  per-human data  ———
GVI::GVI (const Transmission::TransmissionModel& tm) :
    initialInsecticide( 0.0 )   // start with no insecticide (for monitoring)
{
    //TODO: we shouldn't really have vector intervention data (this class) if there's no vector
    // model, should we? Allocate dynamically or based on model?
    const Transmission::VectorModel* vt = dynamic_cast<const Transmission::VectorModel*>(&tm);
    if( vt != 0 ){
        const GVIParams& params = vt->getGVIParams();
        if( params.decay.get() == 0 )
            return;     // no intervention
        // Varience factor of decay is sampled once per human: human is assumed
        // to account for most variance.
        decayHet = params.decay->hetSample();
    }
}

void GVI::deploy(const GVIParams& params) {
    deployTime = TimeStep::simulation;
}

double GVI::relativeAttractiveness(const GVIAnophelesParams& params) const{
    double effect = (1.0 - params._relativeAttractiveness *
            getEffectSurvival(*params.base));
    return params.byProtection( effect );
}

double GVI::preprandialSurvivalFactor(const GVIAnophelesParams& params) const{
    double effect = (1.0 - params._preprandialKillingEffect *
            getEffectSurvival(*params.base));
    return params.byProtection( effect );
}

double GVI::postprandialSurvivalFactor(const GVIAnophelesParams& params) const{
    double effect = (1.0 - params._postprandialKillingEffect *
            getEffectSurvival(*params.base));
    return params.byProtection( effect );
}

} }
