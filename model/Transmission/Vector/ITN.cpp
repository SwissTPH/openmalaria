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
#include "R_nmath/qnorm.h"

namespace OM { namespace Transmission {
    using util::random::poisson;

double ITNParams::init( const scnXml::ITNDescription& elt) {
    initialInsecticide.setParams( elt.getInitialInsecticide() );
    const double maxProp = 0.999;       //NOTE: this could be exposed in XML, but probably doesn't need to be
    maxInsecticide = R::qnorm5(maxProp, initialInsecticide.getMu(), initialInsecticide.getSigma(), true, false);
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
    const ITNParams& params,
    const scnXml::ITNDescription::AnophelesParamsType& elt,
    double proportionUse)
{
    _relativeAttractiveness.init( params, elt.getDeterrency() );
    _preprandialKillingEffect.init( params, elt.getPreprandialKillingEffect(), false );
    _postprandialKillingEffect.init( params, elt.getPostprandialKillingEffect(), true );
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
void ITNAnophelesParams::RelativeAttractiveness::init(const ITNParams& params, const scnXml::ITNDeterrency& elt){
    double HF = elt.getHoleFactor();
    double PF = elt.getInsecticideFactor();
    double IF = elt.getInteractionFactor();
    holeScaling = elt.getHoleScalingFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    if( !(holeScaling>=0.0 && insecticideScaling>=0.0) ){
        throw util::xml_scenario_error("ITN.description.anophelesParams.deterrency: expected scaling factors to be non-negative");
    }
    
    /* We need to ensure the relative availability is non-negative. However,
     * since it's an exponentiated value, it always will be.
     * 
     * If don't want nets to be able to increase transmission, the following
     * limits could also be applied. In general, however, there is no reason
     * nets couldn't make individuals more attractive to mosquitoes.
     * 
     * To ensure relative availability is at most one: relative availability is
     *  exp( log(HF)*h + log(PF)*p + log(IF)*h*p )
     * where HF, PF and IF are the hole, insecticide and interaction factors
     * respectively, with h and p defined as:
     *  h=exp(-holeIndex*holeScalingFactor),
     *  p=1−exp(-insecticideContent*insecticideScalingFactor).
     * We therefore need to ensure that:
     *  log(HF)*h + log(PF)*p + log(IF)*h*p ≤ 0
     * 
     * As with the argument below concerning limits of the killing effect
     * parameters, h and p will always be in the range [0,1] and p ≤ pmax.
     * We can then derive some bounds for HF and PF:
     *  log(HF) ≤ 0
     *  log(PF)×pmax = log(PF^pmax) ≤ 0
     *  log(HF) + (log(PF)+log(IF))×pmax = log(HF×(PF×IF)^pmax) ≤ 0
     * or equivalently
     *  HF ∈ (0,1]
     *  PF^pmax ∈ (0,1]
     *  HF×(PF×IF)^pmax ∈ (0,1]
     *
     * Weaker limits would not be sufficient, as with the argument for the
     * limits of killing effect arguments below. */
#ifdef WITHOUT_BOINC
    // Print out a warning if nets may increase transmission, but only in
    // non-BOINC mode, since it is not unreasonable and volunteers often
    // mistake this kind of warning as indicating a problem.
    double pmax = 1.0-exp(-params.maxInsecticide*insecticideScaling);
    if( !( HF > 0.0 && PF > 0.0 && IF > 0.0 &&
            HF <= 1.0 && pow(PF,pmax) <= 1.0 && HF*pow(PF*IF,pmax) <= 1.0 ) )
    {
        cerr << "Note: since the following bounds are not met, the ITN may make humans more\n";
        cerr << "attractive to mosquitoes than they would be without a net.\n";
        cerr << "This note is only shown by non-BOINC executables.\n";
        cerr << "ITN.description.anophelesParams.deterrency: bounds not met:\n";
        if( !(HF>0.0) )
            cerr << "  holeFactor>0\n";
        if( !(PF>0.0) )
            cerr << "  insecticideFactor>0\n";
        if( !(IF>0.0) )
            cerr << "  interactionFactor>0\n";
        if( !(HF<=1.0) )
            cerr << "  holeFactor≤1\n";
        if( !(pow(PF,pmax)<=1.0) )
            cerr << "  insecticideFactor^"<<pmax<<"≤1\n";
        if( !(HF*pow(PF*IF,pmax)<=1.0) )
            cerr << "  holeFactor×(insecticideFactor×interactionFactor)^"<<pmax<<"≤1\n";
        cerr.flush();
    }
#endif
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
void ITNAnophelesParams::SurvivalFactor::init(const ITNParams& params, const scnXml::ITNKillingEffect& elt, bool postPrandial){
    BF = elt.getBaseFactor();
    HF = elt.getHoleFactor();
    PF = elt.getInsecticideFactor();
    IF = elt.getInteractionFactor();
    holeScaling = elt.getHoleScalingFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    invBaseSurvival = 1.0 / (1.0 - BF);
    if( !( BF >= 0.0 && BF < 1.0) ){
        ostringstream msg;
        msg << "ITN.description.anophelesParams." << (postPrandial?"post":"pre") << "killingFactor: expected baseFactor to be in range [0,1]";
        throw util::xml_scenario_error( msg.str() );
    }
    if( !(holeScaling>=0.0 && insecticideScaling>=0.0) ){
        ostringstream msg;
        msg << "ITN.description.anophelesParams." << (postPrandial?"post":"pre") << "killingFactor: expected scaling factors to be non-negative";
        throw util::xml_scenario_error( msg.str() );
    }
    
    /* We want the calculated survival factor (1−K)/(1−BF) to be in the range
    [0,1] where K is the killing factor: K=BF+HF×h+PF×p+IF×h×p, with h and p
    defined as: h=exp(-holeIndex×holeScalingFactor),
    p=1−exp(-insecticideContent×insecticideScalingFactor). 
    
    By their nature, holeIndex ≥ 0 and insecticideContent ≥ 0. We restrict:
        holeScalingFactor ≥ 0
        insecticideScalingFactor ≥ 0
    Which implies both h and p lie in the range [0,1]. We also know the base
    survival factor, 1−BF, is in the range [0,1].
    
    To make sure the survival factor is not negative we need (1−K)/(1−BF) ≥ 0.
    Since 1−BF > 0 we need 1−K ≥ 0, which, substituting K, gives us
        BF + HF×h + PF×p + IF×h×p ≤ 1	(1)
    We also want to make sure the survival factor is not greater than one (since
    nets shouldn't increase mosquito survival), (1−K)/(1−BF) ≤ 1 or equivalently
    1-K ≤ 1-BF or K ≥ BF, which, substituting K, yields
        HF×h + PF×p + IF×h×p ≥ 0		(2)
    
    Lets derive some limits on HF, PF and IF such that the above inequalities
    (1) and (2) are satisfied.
    
    A net can theoretically be unholed (holeIndex=0 ⇒ h=1) and have no
    insecticide (thus have p=0). Substituting these values in (1) and (2) yields:
        BF + HF ≤ 1	(3)
        HF ≥ 0		(4)
    
    The maximum value for p depends on the maximum insecticide content; denote
    pmax = max(p). Note that holeIndex has no finite maximum; thus, although
    for any finite value of holeIndex, h > 0, there is no h₀ > 0 s.t. for all
    values of holeIndex h ≥ h₀. For the limiting case of a tattered but
    insecticide-saturated net our parameters are therefore p=pmax, h=0:
        BF + PF×pmax ≤ 1	(5)
        PF×pmax ≥ 0			(6)
    (Assuming pmax > 0, (6) is equivalent to PF ≥ 0.)
    
    Consider a net saturated with insecticide (p=pmax) and without holes (h=1):
        BF + HF + (PF+IF)×pmax ≤ 1	(7)
        HF + (PF+IF)×pmax ≥ 0		(8)
    
    The opposite extreme (the limiting case of a decayed net with no remaining
    insecticide and a large number of holes) yields only BF ≤ 1 which we already
    know.
    
    Some of the above examples of nets may be unlikely, but there is only one
    restriction in our model making any of these cases impossible: some
    insecticide must have been lost by the time any holes occur. We ignore this
    since its effect is likely small, and thus all of the above are required to
    keep the survival factor in the range [0,1]. Further, these six inequalities
    (3) - (8) are sufficient to keep the survival factor within [0,1] since h
    and p are non-negative and act linearly in (1) and (2).
    
    From the definition of p, we always have pmax ≤ 1, so substituting pmax=1
    in (5) - (8) gives us bounds which imply our requirement, however if pmax
    is finite they are stricter than necessary. Since insecticideScalingFactor
    is constant, max(p) coincides with max(insecticideContent) which, since
    insecticide content only decays over time, coincides with the maximum
    initial insecticide content, Pmax. Since the initial insecticide content is
    sampled from a normal distribution in our model it should have no finite
    maximum, thus implying we cannot achieve more relaxed bounds than (5) - (8)
    when pmax=1 (unless the standard deviation of our normal distribution is 1).
    We would however like to impose less strict bounds than these, thus we
    impose a maximum value on the initial insecticide content, Pmax, such that
    the probability of sampling a value from our parameterise normal
    distribution greater than Pmax is 0.001. */
    double pmax = 1.0-exp(-params.maxInsecticide*insecticideScaling);
    if( !( BF+HF <= 1.0 && HF >= 0.0
        && BF+PF*pmax <= 1.0 && PF*pmax >= 0.0
        && BF+HF+(PF+IF)*pmax <= 1.0 && HF+(PF+IF)*pmax >= 0.0 ) )
    {
        ostringstream msg;
        msg << "ITN.description.anophelesParams." << (postPrandial?"post":"pre") << "killingFactor: bounds not met:";
        if( !(BF+HF<=1.0) )
            msg << " baseFactor+holeFactor≤1";
        if( !(HF>=0.0) )
            msg << " holeFactor≥0";
        if( !(BF+PF*pmax<=1.0) )
            msg << " baseFactor+"<<pmax<<"×insecticideFactor≤1";
        if( !(PF*pmax>=0.0) )
            msg << " insecticideFactor≥0";      // if this fails, we know pmax>0 (since it is in any case non-negative) — well, or an NaN
        if( !(PF+HF+(PF+IF)*pmax<=1.0) )
            msg << " baseFactor+holeFactor+"<<pmax<<"×(insecticideFactor+interactionFactor)≤1";
        if( !(HF+(PF+IF)*pmax>=0.0) )
            msg << " holeFactor+"<<pmax<<"×(insecticideFactor+interactionFactor)≥0";
        throw util::xml_scenario_error( msg.str() );
    }
}
double ITNAnophelesParams::RelativeAttractiveness::relativeAttractiveness( double holeIndex, double insecticideContent )const {
    double holeComponent = exp(-holeIndex*holeScaling);
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double relAvail = exp( lHF*holeComponent + lPF*insecticideComponent + lIF*holeComponent*insecticideComponent );
    assert( relAvail>=0.0 );
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

void ITN::deploy(TimeStep now, const ITNParams& params) {
    deployTime = now;
    disposalTime = now + params.attritionOfNets->sampleAgeOfDecay();
    nHoles = 0;
    holeIndex = 0.0;
    // this is sampled independently: initial insecticide content doesn't depend on handling
    initialInsecticide = params.initialInsecticide.sample();
    if( initialInsecticide < 0.0 )
        initialInsecticide = 0.0;	// avoid negative samples
    if( initialInsecticide > params.maxInsecticide )
        initialInsecticide = params.maxInsecticide;
    
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
