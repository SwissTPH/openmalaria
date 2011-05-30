/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Transmission/Vector/ITN.h"
#include "util/random.h"
#include "util/errors.h"

namespace OM { namespace Transmission {
    using util::random::poisson;

double ITNParams::init( const scnXml::ITNDescription& elt) {
    initialInsecticide.setParams( elt.getInitialInsecticide() );
    holeRate.setParams( elt.getHoleRate() );
    ripRate.setParams( elt.getRipRate() );
    ripFactor = elt.getRipFactor().getValue();
    insecticideDecay = DecayFunction::makeObject( elt.getInsecticideDecay(), "ITNDescription.insecticideDecay" );
    double propUse = elt.getUsage().getValue();
    if( !( propUse >= 0.0 && propUse <= 1.0 ) ){
        throw util::xml_scenario_error("ITN.description.proportionUse: must be within range [0,1]");
    }
    return propUse;
}

void ITNAnophelesParams::init(
    const scnXml::ITNDescription::AnophelesParamsType& elt, double proportionUse)
{
    _relativeAvailability.init( elt.getRelativeAvailability() );
    _preprandialKillingEffect.init( elt.getPreprandialKillingEffect() );
    _postprandialKillingEffect.init( elt.getPostprandialKillingEffect() );
    // Nets only affect people while they're indoors and using the net.
    // TODO: This could perhaps be extracted from ITN code since it also
    // affects IRS and mosquito deterrents.
    assert( proportionUse >= 0.0 && proportionUse <= 1.0 );
    double propBitingIndoors = elt.getProportionBitingIndoors().getValue();
    if( !(propBitingIndoors >= 0.0 && propBitingIndoors <= 1.0 ) ){
        throw util::xml_scenario_error(
            "ITN.description.anophelesParams.proportionBitingIndoors: must be within range [0,1]"
        );
    }
    proportionProtected = proportionUse * propBitingIndoors;
    proportionUnprotected = 1.0 - proportionProtected;
}

inline bool inRange01( double x ){
    return x>=0.0 && x<= 1.0;
}
void ITNAnophelesParams::RelativeAvailability::init(const scnXml::ITNAvailEffect& elt){
    // Factors are adjusted such that we don't need to take 1 minus
    // the exponential of neg. scaled insecticide content.
    basePlusInsecticide = elt.getInsecticideFactor();
    holePlusInteraction = elt.getHoleFactor() + elt.getInteractionFactor();
    negInsecticide = -elt.getInsecticideFactor();
    negInteraction = -elt.getInteractionFactor();
    holeScaling = elt.getHoleScalingFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    if( !(holeScaling>=0.0 && insecticideScaling>=0.0) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams: expected scaling factors to be non-negative");
    }
    
    // We want the calculated effect to be in the range [0,1].
    // Since holeIndex*holeScaling and insecticideContent*insecticideScaling are
    // both non-negative, both exponentials are in the range (0,1]. Considering
    // extreme values for the exponentials and realising that their product is
    // no greater than either alone leads us to the following requirements:
    // hole and insecticide factors are in the range [0,1], and
    // hole+insecticide+interaction factors are in the range [0,1].
    if( !(inRange01(elt.getInsecticideFactor()) &&
            inRange01(elt.getHoleFactor()) &&
            inRange01(holePlusInteraction+elt.getInsecticideFactor())) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams: expected hole factor, insecticide factor and hole+insecticide+interaction factors all to be in range [0,1]");
    }
}
void ITNAnophelesParams::KillingEffect::init(const scnXml::ITNKillingEffect& elt){
    ITNAnophelesParams::RelativeAvailability::init( elt );
    // We want to assert that (1-killingFactor) <= (1-baseFactor). This turns
    // out to be the same as requiring the effect >= 0 as above.
    
    if( !(elt.getBaseFactor() >= 0.0 && elt.getBaseFactor() <= 1.0) ){
        // is a probability so must be in range [0,1]
        throw util::xml_scenario_error("ITN.description.anophelesParams: expected baseFactor to be in range [0,1]");
    }
    basePlusInsecticide += elt.getBaseFactor();
    invBaseSurvival = 1.0 / (1.0 - elt.getBaseFactor());
}
double ITNAnophelesParams::RelativeAvailability::effect( double holeIndex, double insecticideContent )const {
    double holeComponent = exp(-holeIndex*holeScaling);
    // The model formula uses 1 minus the exponential here. We adjust
    // the factors in the constructor to achieve the same effect.
    double insecticideComponent = exp(-insecticideContent*insecticideScaling);
    double unadjusted = basePlusInsecticide
        + holePlusInteraction*holeComponent
        + negInsecticide*insecticideComponent
        + negInteraction*holeComponent*insecticideComponent;
    assert( inRange01(unadjusted) );
    return unadjusted;
}
double ITNAnophelesParams::KillingEffect::effect( double holeIndex, double insecticideContent )const {
    double killingEffect = RelativeAvailability::effect( holeIndex, insecticideContent );
    double survivalFactor = (1.0 - killingEffect) * invBaseSurvival;
    assert( inRange01(survivalFactor) );
    return survivalFactor;
}

void ITN::deploy(const ITNParams& params) {
    deployTime = TimeStep::simulation;
    nHoles = 0;
    holeIndex = 0.0;
    initialInsecticide = params.initialInsecticide.sample();
    
    holeRate = params.holeRate.sample();
    ripRate = params.ripRate.sample();
    insecticideDecayHet = params.insecticideDecay->hetSample();
}

void ITN::update(const ITNParams& params){
    if( deployTime != TimeStep::never ){
        int newHoles = poisson( holeRate );
        nHoles += newHoles;
        holeIndex += newHoles + params.ripFactor * poisson( nHoles * ripRate );
    }
}

double ITN::relativeAvailability(const ITNAnophelesParams& params) const{
    return params.relativeAvailability( holeIndex, initialInsecticide * params.base->insecticideDecay->eval (TimeStep::simulation - deployTime, insecticideDecayHet));
}

double ITN::preprandialSurvivalFactor(const ITNAnophelesParams& params) const{
    return params.preprandialSurvivalFactor( holeIndex, initialInsecticide * params.base->insecticideDecay->eval (TimeStep::simulation - deployTime, insecticideDecayHet));
}

double ITN::postprandialSurvivalFactor(const ITNAnophelesParams& params) const{
    return params.postprandialSurvivalFactor( holeIndex, initialInsecticide * params.base->insecticideDecay->eval (TimeStep::simulation - deployTime, insecticideDecayHet));
}

} }
