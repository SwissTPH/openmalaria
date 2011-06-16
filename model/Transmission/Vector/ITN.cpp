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
    attritionOfNets = DecayFunction::makeObject( elt.getAttritionOfNets(), "ITNDescription.attritionOfNets" );
    double propUse = elt.getUsage().getValue();
    if( !( propUse >= 0.0 && propUse <= 1.0 ) ){
        throw util::xml_scenario_error("ITN.description.proportionUse: must be within range [0,1]");
    }
    return propUse;
}

void ITNAnophelesParams::init(
    const scnXml::ITNDescription::AnophelesParamsType& elt, double proportionUse)
{
    _relativeAttractiveness.init( elt.getDeterrency() );
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
ITNAnophelesParams::RelativeAttractiveness::RelativeAttractiveness() :
    lHF( numeric_limits< double >::signaling_NaN() ),
    lPF( numeric_limits< double >::signaling_NaN() ),
    lIF( numeric_limits< double >::signaling_NaN() ),
    holeScaling( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() )
{}
void ITNAnophelesParams::RelativeAttractiveness::init(const scnXml::ITNDeterrency& elt){
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
    
    h and p will always be in the range [0,1], so the following limits are
    sufficient to keep the relative availability in the range [0,1]:
    log(HF) ≤ 0, log(PF) ≤ 0 and log(HF)+log(PF)+log(IF) = log(HF×PF×IF) ≤ 0,
    or equivalently HF∈(0,1], PF∈(0,1] and HF×PF×IF∈(0,1].
    Weaker limits would not be sufficient, as with the argument for the limits
    of killing effect arguments below.
    */
    if( !( HF > 0.0 && HF <= 1.0 && PF > 0.0 && PF <= 1.0 &&
            HF*PF*IF > 0.0 && HF*PF*IF <= 1.0 ) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams.deterrency: "
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
    [0,1] where K is the killing factor: K=BF+HF×h+PF×p+IF×h×p, with h and p
    defined as: h=exp(-holeIndex×holeScalingFactor),
    p=1−exp(-insecticideContent×insecticideScalingFactor).
    
    Since 1−BF > 0 we need 1−K ≥ 0 or equivalently BF+HF×h+PF×p+IF×h×p ≤ 1.
    We also need 1-K ≤ 1-BF, or equivalently K ≥ BF which reduces to
    HF×h + PF×p + IF×h×p ≥ 0.
    
    Initially h is 1 (when holeIndex=0), and since holeIndex increases towards
    infinity with net age, h will tend to 0 (assuming the scaling factor is
    positive). p will tend towards 0 as the insecticide content tends towards
    0, however it's initial value p₀ depends on the initial content and the
    scaling factor. We know 0 ≤ p ≤ p₀ ≤ 1 but not the exact value of p₀ since
    insecticideContent is sampled from a normal distribution per net; in theory
    though the maximum sample is unlimited and thus p₀ can take any value less
    than 1.
    
    Substituting the initial values for h and p we require:
    BF + HF + (PF+IF)×p₀ ≤ 1, and HF + (PF+IF)×p₀ ≥ 0.
    Since p₀ may be as low as 0, we require HF ≥ 0 and BF + HF ≤ 1. The extreme
    when both h and p tend to zero doesn't tell us much, but it is also possible
    for h and p to take values in between. Exactly what limit would be needed
    in all these cases isn't obvious, but certainly requiring PF ≥ 0 and
    BF + PF×p₀ ≤ 1 (which assumes h→0 while p=p₀) is sufficient.
    
    Since we want limits without p₀, observe that the initial insecticide
    content is sampled from a normal distribution and thus has no maximum
    bound, implying that p₀ can take any value less than 1, so substituting p₀
    with the value 1 is the only way to guarantee that (1−K)/(1−BF) ∈ [0,1].
    */
    if( !( BF+HF <= 1.0 && BF+PF <= 1.0 && BF+HF+PF+IF <= 1.0 &&
        HF >= 0.0 && PF >= 0.0 && HF+PF+IF >= 0.0 ) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams.*killingFactor: "
        "bounds not met: BF+HF ≤ 1, BF+PF ≤ 1, BF+HF+PF+IF ≤ 1, HF ≥ 0, PF ≥ 0 and HF+PF+IF ≥ 0" );
    }
}
double ITNAnophelesParams::RelativeAttractiveness::relativeAttractiveness( double holeIndex, double insecticideContent )const {
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
    disposalTime = TimeStep::simulation + params.attritionOfNets->sampleAgeOfDecay();
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
        if( TimeStep::simulation >= disposalTime ){
            deployTime = TimeStep::never;
        }
        int newHoles = poisson( holeRate );
        nHoles += newHoles;
        holeIndex += newHoles + params.ripFactor * poisson( nHoles * ripRate );
    }
}

double ITN::relativeAttractiveness(const ITNAnophelesParams& params) const{
    return params.relativeAttractiveness( holeIndex, getInsecticideContent(*params.base) );
}

double ITN::preprandialSurvivalFactor(const ITNAnophelesParams& params) const{
    return params.preprandialSurvivalFactor( holeIndex, getInsecticideContent(*params.base) );
}

double ITN::postprandialSurvivalFactor(const ITNAnophelesParams& params) const{
    return params.postprandialSurvivalFactor( holeIndex, getInsecticideContent(*params.base) );
}

} }
