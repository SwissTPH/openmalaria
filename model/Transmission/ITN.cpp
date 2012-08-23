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

#include "Transmission/ITN.h"
//TODO: we shouldn't have a dependency on the vector/transmission model class
//here; currently it's a work-around for ITN parameters not always being present.
#include "Transmission/VectorModel.h"
#include "util/random.h"
#include "util/errors.h"
#include "R_nmath/qnorm.h"
#include <cmath>

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
    if (elt.getDeterrency().present())
        _relativeAttractiveness = shared_ptr<RelativeAttractiveness>(new RADeterrency( params, elt.getDeterrency().get() ));
    else{
        assert (elt.getTwoStageDeterrency().present());
        _relativeAttractiveness = shared_ptr<RelativeAttractiveness>(new RATwoStageDeterrency( params, elt.getTwoStageDeterrency().get() ));
    }
    _preprandialKillingEffect.init( params, elt.getPreprandialKillingEffect(), "ITN.description.anophelesParams.preprandialKillingFactor" );
    _postprandialKillingEffect.init( params, elt.getPostprandialKillingEffect(), "ITN.description.anophelesParams.postprandialKillingFactor" );
    // Nets only affect people while they're using the net. NOTE: we may want
    // to revise this at some point (heterogeneity, seasonal usage patterns).
    double propActive = elt.getPropActive();
    assert( proportionUse >= 0.0 && proportionUse <= 1.0 );
    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = proportionUse * propActive;
    proportionUnprotected = 1.0 - proportionProtected;
}

ITNAnophelesParams::RADeterrency::RADeterrency(const ITNParams& params, const scnXml::ITNDeterrency& elt) :
    lHF( numeric_limits< double >::signaling_NaN() ),
    lPF( numeric_limits< double >::signaling_NaN() ),
    lIF( numeric_limits< double >::signaling_NaN() ),
    holeScaling( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() )
{
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
     *  log(PF)×pmax ≤ 0
     *  log(HF) + (log(PF)+log(IF))×pmax = log(HF×(PF×IF)^pmax) ≤ 0
     * or equivalently (assuming pmax>0):
     *  HF ∈ (0,1]
     *  PF ∈ (0,1]
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
            HF <= 1.0 && PF <= 1.0 && HF*pow(PF*IF,pmax) <= 1.0 ) )
    {
        cerr << "Note: since the following bounds are not met, the ITN could make humans more\n";
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
        if( !(PF<=1.0) )
            cerr << "  insecticideFactor≤1\n";
        if( !(HF*pow(PF*IF,pmax)<=1.0) )
            cerr << "  holeFactor×(insecticideFactor×interactionFactor)^"<<pmax<<"≤1\n";
        cerr.flush();
    }
#endif
    lHF = log( HF );
    lPF = log( PF );
    lIF = log( IF );
}
ITNAnophelesParams::RATwoStageDeterrency::RATwoStageDeterrency(
        const OM::Transmission::ITNParams& params,
        const scnXml::TwoStageDeterrency& elt) :
    lPFEntering( numeric_limits< double >::signaling_NaN() ),
    insecticideScalingEntering( numeric_limits< double >::signaling_NaN() )
{
    //TODO: this is just a copy from IRSAnophelesParams::RelativeAttractiveness::init
    // It should be possible to abstract out a lot of this code.
    
    double PF = elt.getEntering().getInsecticideFactor();
    insecticideScalingEntering = elt.getEntering().getInsecticideScalingFactor();
    if( !( PF > 0.0) ){
        ostringstream msg;
        msg << "ITN.description.anophelesParams.twoStageDeterrency.entering: expected insecticideFactor to be positive.";
        //TODO: These constraints were required. But they're too strong.
        // Now need to work out which should still be imposed.
        cerr << msg.str() << endl;
        //throw util::xml_scenario_error( msg.str() );
    }
    
    /* We need to ensure the relative availability is non-negative. However,
     * since it's an exponentiated value, it always will be.
     * 
     * If we don't want ITNs to be able to increase transmission, the following
     * limits could also be applied. In general, however, there is no reason
     * ITNs couldn't make individuals more attractive to mosquitoes.
     * 
     * To ensure relative availability is at most one: relative availability is
     *  exp( log(PF)*p ) = PF^p
     * where PF is the insecticide factor, with p∈[0,1] defined as:
     *  p=1−exp(-insecticideContent*insecticideScalingFactor).
     * We therefore just need PF ≤ 1. */
#ifdef WITHOUT_BOINC
    // Print out a warning if ITNs may increase transmission, but only in
    // non-BOINC mode, since it is not unreasonable and volunteers often
    // mistake this kind of warning as indicating a problem.
    if( !( PF <= 1.0 ) ) {
        cerr << "Note: since the following bounds are not met, the IRS could make humans more\n";
        cerr << "attractive to mosquitoes than they would be without IRS.\n";
        cerr << "This note is only shown by non-BOINC executables.\n";
        cerr << "IRS.description.anophelesParams.deterrency: bounds not met:\n";
        cerr << "  0<insecticideFactor≤1\n";
        cerr.flush();
    }
#endif
    lPFEntering = log( PF );
    
    pAttacking.init( params, elt.getAttacking(), "ITN.description.anophelesParams.twoStageDeterrency.attacking" );
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
void ITNAnophelesParams::SurvivalFactor::init(const OM::Transmission::ITNParams& params, const scnXml::ITNKillingEffect& elt, const char* eltName){
    BF = elt.getBaseFactor();
    HF = elt.getHoleFactor();
    PF = elt.getInsecticideFactor();
    IF = elt.getInteractionFactor();
    holeScaling = elt.getHoleScalingFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    invBaseSurvival = 1.0 / (1.0 - BF);
    if( !( BF >= 0.0 && BF < 1.0) ){
        ostringstream msg;
        msg << eltName << ": expected baseFactor to be in range [0,1]";
        throw util::xml_scenario_error( msg.str() );
    }
    if( !(holeScaling>=0.0 && insecticideScaling>=0.0) ){
        ostringstream msg;
        msg << eltName << ": expected scaling factors to be non-negative";
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
        && BF+PF*pmax <= 1.0 && PF >= 0.0
        && BF+HF+(PF+IF)*pmax <= 1.0 && HF+(PF+IF)*pmax >= 0.0 ) )
    {
        ostringstream msg;
        msg << eltName << ": bounds not met:";
        if( !(BF+HF<=1.0) )
            msg << " baseFactor+holeFactor≤1";
        if( !(HF>=0.0) )
            msg << " holeFactor≥0";
        if( !(BF+PF*pmax<=1.0) )
            msg << " baseFactor+"<<pmax<<"×insecticideFactor≤1";
        if( !(PF>=0.0) )
            msg << " insecticideFactor≥0";      // if this fails, we know pmax>0 (since it is in any case non-negative) — well, or an NaN
        if( !(PF+HF+(PF+IF)*pmax<=1.0) )
            msg << " baseFactor+holeFactor+"<<pmax<<"×(insecticideFactor+interactionFactor)≤1";
        if( !(HF+(PF+IF)*pmax>=0.0) )
            msg << " holeFactor+"<<pmax<<"×(insecticideFactor+interactionFactor)≥0";
        //TODO: These constraints were required. But they're too strong.
        // Now need to work out which should still be imposed.
        cerr << msg.str() << endl;
        //throw util::xml_scenario_error( msg.str() );
    }
}
double ITNAnophelesParams::RADeterrency::relativeAttractiveness( double holeIndex, double insecticideContent )const {
    double holeComponent = exp(-holeIndex*holeScaling);
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double relAvail = exp( lHF*holeComponent + lPF*insecticideComponent + lIF*holeComponent*insecticideComponent );
    //TODO: limits
    //assert( relAvail>=0.0 );
    if (relAvail < 0)
        relAvail = 0.0;
    return relAvail;
}
double ITNAnophelesParams::RATwoStageDeterrency::relativeAttractiveness(
        double holeIndex, double insecticideContent )const
{
    // This is essentially a combination of the relative attractiveness as used
    // by IRS and a killing factor.
    
    // Note that an alternative, simpler, model could have been used, but was
    // not for consistency with other models. Alternative (here we don't take
    // the logarithm of PF):
    // pEnt = 1 - PFEntering × insecticideComponent
    
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScalingEntering);
    double pEnt = exp( lPFEntering*insecticideComponent );
    assert( pEnt >= 0.0 );
    
    double rel_pAtt = pAttacking.rel_pAtt( holeIndex, insecticideContent );
    // normalise: must have 1 when no insecticide and no net (infinite holes):
    //TODO: limits
    if (pEnt * rel_pAtt < 0.0)
        return 0.0;
    return pEnt * rel_pAtt;
}
double ITNAnophelesParams::SurvivalFactor::rel_pAtt( double holeIndex, double insecticideContent )const {
    double holeComponent = exp(-holeIndex*holeScaling);
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double pAtt = BF + HF*holeComponent + PF*insecticideComponent + IF*holeComponent*insecticideComponent;
    //TODO: limits
    //assert( pAtt <= 1.0 );
    return pAtt / BF;
}
double ITNAnophelesParams::SurvivalFactor::survivalFactor( double holeIndex, double insecticideContent )const {
    double holeComponent = exp(-holeIndex*holeScaling);
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double killingEffect = BF + HF*holeComponent + PF*insecticideComponent + IF*holeComponent*insecticideComponent;
    //assert( killingEffect <= 1.0 );
    double survivalFactor = (1.0 - killingEffect) * invBaseSurvival;
    //TODO: limits
    //assert( survivalFactor >= 0.0 );
    if (survivalFactor < 0.0)
        return 0.0;
    else if (survivalFactor > 1.0)
        return 1.0;
    return survivalFactor;
}

ITN::ITN(const TransmissionModel& tm) :
        nHoles( 0 ),
        holeIndex( numeric_limits<double>::signaling_NaN() ),
        initialInsecticide( numeric_limits<double>::signaling_NaN() ),
        holeRate( numeric_limits<double>::signaling_NaN() ),
        ripRate( numeric_limits<double>::signaling_NaN() )
{
    //TODO: we shouldn't really have ITN data (this class) if there's no vector
    // model, should we? Allocate dynamically or based on model?
    const VectorModel* vt = dynamic_cast<const VectorModel*>(&tm);
    if( vt != 0 ){
        const ITNParams& params = vt->getITNParams();
        if( params.insecticideDecay.get() == 0 )
            return;     // no ITNs
        // Net rips and insecticide loss are assumed to co-vary dependent on
        // handling of net. They are sampled once per human: human handling is
        // presumed to be the largest cause of variance.
        util::NormalSample x = util::NormalSample::generate();
        holeRate = params.holeRate.sample(x) * TimeStep::yearsPerInterval;
        ripRate = params.ripRate.sample(x) * TimeStep::yearsPerInterval;
        insecticideDecayHet = params.insecticideDecay->hetSample(x);
    }
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
    if( initialInsecticide > params.maxInsecticide )
        initialInsecticide = params.maxInsecticide;
}

void ITN::update(const ITNParams& params){
    if( deployTime != TimeStep::never ){
        // First use is at age 1, so don't remove until *after* disposalTime to
        // get use over the full duration given by sampleAgeOfDecay().
        if( TimeStep::simulation > disposalTime ){
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
