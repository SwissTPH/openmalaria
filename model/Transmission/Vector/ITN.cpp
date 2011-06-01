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
    // Nets only affect people while they're using the net. NOTE: we may want
    // to revise this at some point (heterogeneity, seasonal usage patterns).
    assert( proportionUse >= 0.0 && proportionUse <= 1.0 );
    proportionProtected = proportionUse;
    proportionUnprotected = 1.0 - proportionProtected;
}

inline bool inRange01( double x ){
    return x>=0.0 && x<= 1.0;
}
ITNAnophelesParams::RelativeAvailability::RelativeAvailability() :
    lHF( numeric_limits< double >::signaling_NaN() ),
    lPF( numeric_limits< double >::signaling_NaN() ),
    lIF( numeric_limits< double >::signaling_NaN() ),
    holeScaling( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() )
{}
void ITNAnophelesParams::RelativeAvailability::init(const scnXml::ITNAvailEffect& elt){
    double HF = elt.getHoleFactor();
    double PF = elt.getInsecticideFactor();
    double IF = elt.getInteractionFactor();
    holeScaling = elt.getHoleScalingFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    if( !(holeScaling>=0.0 && insecticideScaling>=0.0) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams: expected scaling factors to be non-negative");
    }
    
    /* We need to ensure the relative availability is in the range [0,1]. It's
    an exponentiated value, so we require log(HF)*h + log(PF)*p + log(IF)*h*p ≤ 0
    where HF, PF and IF are the hole, insecticide and interaction factors
    respectively, with h and p defined as: h=exp(-holeIndex*holeScalingFactor),
    p=1−exp(-insecticideContent*insecticideScalingFactor).
    
    h and p will always be in the range [0,1], so we need log(HF) ≤ 0,
    log(PF) ≤ 0 and log(HF)+log(PF)+log(IF) = log(HF×PF×IF) ≤ 0.
    This implies that we need HF∈(0,1], PF∈(0,1] and HF×PF×IF∈(0,1].
    */
    if( !( HF > 0.0 && HF <= 1.0 && PF > 0.0 && PF <= 1.0 &&
            HF*PF*IF > 0.0 && HF*PF*IF <= 1.0 ) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams.relativeAvailability: "
        "bounds not met: HF∈(0,1], PF∈(0,1] and HF×PF×IF∈(0,1]" );
    }
    lHF = log( HF );
    lPF = log( PF );
    lIF = log( IF );
}
ITNAnophelesParams::SurvivalFactor::SurvivalFactor() :
    BF( numeric_limits< double >::signaling_NaN() ),
    HF( numeric_limits< double >::signaling_NaN() ),
    PF( numeric_limits< double >::signaling_NaN() ),
    IF( numeric_limits< double >::signaling_NaN() ),
    holeScaling( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() ),
    invBaseSurvival( numeric_limits< double >::signaling_NaN() )
{}
void ITNAnophelesParams::SurvivalFactor::init(const scnXml::ITNKillingEffect& elt){
    BF = elt.getBaseFactor();
    HF = elt.getHoleFactor();
    PF = elt.getInsecticideFactor();
    IF = elt.getInteractionFactor();
    holeScaling = elt.getHoleScalingFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    invBaseSurvival = 1.0 / (1.0 - BF);
    if( !( BF >= 0.0 && BF < 1.0) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams: expected baseFactor to be in range [0,1]");
    }
    if( !(holeScaling>=0.0 && insecticideScaling>=0.0) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams: expected scaling factors to be non-negative");
    }
    
    /* We want the calculated survival factor (1−K)/(1−BF) to be in the range
    [0,1] where K is the killing factor: K=BF+HF*h+PF*p+IF*h*p, with h and p
    defined as: h=exp(-holeIndex*holeScalingFactor),
    p=1−exp(-insecticideContent*insecticideScalingFactor).
    
    Since 1−BF > 0 we need 1−K ≥ 0 or equivalently BF+HF*h+PF*p+IF*h*p ≤ 1.
    Since holeIndex*holeScaling and insecticideContent*insecticideScaling are
    both non-negative, both h and p are in the range (0,1]. Considering
    extreme values for the exponentials and realising that their product is
    no greater than either alone leads us to the following requirements:
    BF+HF ≤ 1, BF+PF ≤ 1 and BF+HF+PF+IF ≤ 1.
    
    We also need 1-K ≤ 1-BF, or equivalently K ≥ BF which reduces to
    HF*h + PF*p + IF*h*p ≥ 0. Similar to the above, this gives us the requirement
    that HF ≥ 0, PF ≥ 0 and HF+PF+IF ≥ 0.
    */
    if( !( BF+HF <= 1.0 && BF+PF <= 1.0 && BF+HF+PF+IF <= 1.0 &&
        HF >= 0.0 && PF >= 0.0 && HF+PF+IF >= 0.0 ) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams.*killingFactor: "
        "bounds not met: BF+HF ≤ 1, BF+PF ≤ 1, BF+HF+PF+IF ≤ 1, HF ≥ 0, PF ≥ 0 and HF+PF+IF ≥ 0" );
    }
}
double ITNAnophelesParams::RelativeAvailability::relativeAvailability( double holeIndex, double insecticideContent )const {
    double holeComponent = exp(-holeIndex*holeScaling);
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double relAvail = exp( lHF*holeComponent + lPF*insecticideComponent + lIF*holeComponent*insecticideComponent );
    assert( inRange01(relAvail) );
    return relAvail;
}
double ITNAnophelesParams::SurvivalFactor::survivalFactor( double holeIndex, double insecticideContent )const {
    double holeComponent = exp(-holeIndex*holeScaling);
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double killingEffect = BF + HF*holeComponent + PF*insecticideComponent + IF*holeComponent*insecticideComponent;
    assert( inRange01(killingEffect) );
    double survivalFactor = (1.0 - killingEffect) * invBaseSurvival;
    assert( inRange01(survivalFactor) );
    return survivalFactor;
}

void ITN::deploy(const ITNParams& params) {
    deployTime = TimeStep::simulation;
    nHoles = 0;
    holeIndex = 0.0;
    // this is sampled independently: initial insecticide content doesn't depend on handling
    initialInsecticide = params.initialInsecticide.sample();
    if( initialInsecticide < 0.0 )
        initialInsecticide = 0.0;	// avoid negative samples
    
    // net rips and insecticide loss are assumed to co-vary dependent on handling of net
    util::NormalSample x = util::NormalSample::generate();
    holeRate = params.holeRate.sample(x) * TimeStep::yearsPerInterval;
    ripRate = params.ripRate.sample(x) * TimeStep::yearsPerInterval;
    insecticideDecayHet = params.insecticideDecay->hetSample(x);
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
