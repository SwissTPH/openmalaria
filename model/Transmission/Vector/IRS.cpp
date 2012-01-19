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

#include "Transmission/Vector/IRS.h"
//TODO: we shouldn't have a dependency on the vector/transmission model class
//here; currently it's a work-around for IRS parameters not always being present.
#include "Transmission/Vector/VectorTransmission.h"
#include "util/random.h"
#include "util/errors.h"
#include "R_nmath/qnorm.h"
#include <cmath>

namespace OM { namespace Transmission {
    using util::random::poisson;

void IRSParams::init( const scnXml::IRSDescription& elt) {
    simpleModel = false;
    initialInsecticide.setParams( elt.getInitialInsecticide() );
    const double maxProp = 0.999;       //NOTE: this could be exposed in XML, but probably doesn't need to be
    maxInsecticide = R::qnorm5(maxProp, initialInsecticide.getMu(), initialInsecticide.getSigma(), true, false);
    insecticideDecay = DecayFunction::makeObject( elt.getInsecticideDecay(), "IRS.description.insecticideDecay" );
}
void IRSParams::init( const scnXml::IRSSimpleDescription& elt) {
    simpleModel = true;
    insecticideDecay = DecayFunction::makeObject( elt.getDecay(), "IRS.simpleDescription.decay" );
}

void IRSAnophelesParams::init(
    const IRSParams& params,
    const scnXml::IRSDescription::AnophelesParamsType& elt)
{
    assert( !params.simpleModel );
    _relativeAttractiveness.init( params, elt.getDeterrency() );
    _preprandialKillingEffect.init( params, elt.getPreprandialKillingEffect(), false );
    _postprandialKillingEffect.init( params, elt.getPostprandialKillingEffect(), true );
    // Simpler version of ITN usage/action:
    double propActive = elt.getPropActive();
    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = propActive;
    proportionUnprotected = 1.0 - proportionProtected;
}
void IRSAnophelesParams::init(
    const IRSParams& params,
    const scnXml::IRSSimpleDescription::AnophelesParamsType& elt)
{
    assert( params.simpleModel );
    _relativeAttractiveness.oldDeterrency( elt.getDeterrency().getValue() );
    _preprandialKillingEffect.oldEffect( elt.getPreprandialKillingEffect().getValue() );
    _postprandialKillingEffect.oldEffect( elt.getPostprandialKillingEffect().getValue() );
    // Simpler version of ITN usage/action:
    double propActive = elt.getPropActive();
    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = propActive;
    proportionUnprotected = 1.0 - proportionProtected;
}

inline bool inRange01( double x ){
    return x>=0.0 && x<= 1.0;
}
IRSAnophelesParams::RelativeAttractiveness::RelativeAttractiveness() :
    lPF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() )
{}
void IRSAnophelesParams::RelativeAttractiveness::init(const IRSParams& params, const scnXml::IRSDeterrency& elt){
    double PF = elt.getInsecticideFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    if( !( PF > 0.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams.relativeAvailability: expected insecticideFactor to be positive.";
        throw util::xml_scenario_error( msg.str() );
    }
    
    /* We need to ensure the relative availability is non-negative. However,
     * since it's an exponentiated value, it always will be.
     * 
     * If we don't want IRS to be able to increase transmission, the following
     * limits could also be applied. In general, however, there is no reason
     * IRS couldn't make individuals more attractive to mosquitoes.
     * 
     * To ensure relative availability is at most one: relative availability is
     *  exp( log(PF)*p ) = PF^p
     * where PF is the insecticide factor, with p∈[0,1] defined as:
     *  p=1−exp(-insecticideContent*insecticideScalingFactor).
     * We therefore just need PF ≤ 1. */
#ifdef WITHOUT_BOINC
    // Print out a warning if IRS may increase transmission, but only in
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
    lPF = log( PF );
}
IRSAnophelesParams::SurvivalFactor::SurvivalFactor() :
    BF( numeric_limits< double >::signaling_NaN() ),
    PF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() ),
    invBaseSurvival( numeric_limits< double >::signaling_NaN() )
{}
void IRSAnophelesParams::SurvivalFactor::init(const IRSParams& params, const scnXml::IRSKillingEffect& elt, bool postPrandial){
    BF = elt.getBaseFactor();
    PF = elt.getInsecticideFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    invBaseSurvival = 1.0 / (1.0 - BF);
    if( !( BF >= 0.0 && BF < 1.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected baseFactor to be in range [0,1]";
        throw util::xml_scenario_error( msg.str() );
    }
    if( !(insecticideScaling>=0.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected scaling factor to be non-negative";
        throw util::xml_scenario_error( msg.str() );
    }
    
    /* We want the calculated survival factor (1−K)/(1−BF) to be in the range
    [0,1] where K is the killing factor: K=BF+PF×p, with p defined as:
    p=1−exp(-insecticideContent×insecticideScalingFactor).
    
    By its nature, insecticideContent ≥ 0. We restrict
        insecticideScalingFactor ≥ 0
    which implies that p lies in the range [0,1]. We also know the base
    survival factor, 1−BF, is in the range [0,1].
    
    To make sure the survival factor is not negative we need (1−K)/(1−BF) ≥ 0.
    Since 1−BF > 0 we need 1−K ≥ 0, which, substituting K, gives us
        BF + PF×p ≤ 1	(1)
    We also want to make sure the survival factor is not greater than one (since
    IRS shouldn't increase mosquito survival), (1−K)/(1−BF) ≤ 1 or equivalently
    1-K ≤ 1-BF or K ≥ BF, which, substituting K, yields
        PF×p ≥ 0		(2)
    Since p ≥ 0, ensuring PF is non-negative is enough for (2).
    
    The maximum value for p depends on the maximum insecticide content; denote
    pmax = max(p). In this extreme case (1) becomes:
        BF + PF×pmax ≤ 1	(3)
    
    As with the ITN model, we impose a maximum value on the initial insecticide
    content, Pmax, such that the probability of sampling a value from our
    parameterise normal distribution greater than Pmax is 0.001. */
    double pmax = 1.0-exp(-params.maxInsecticide*insecticideScaling);
    if( !( PF >= 0.0 && BF+PF*pmax <= 1.0 ) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected insecticideFactor≥0, baseFactor+"
            <<pmax<<"×insecticideFactor≤1";
        throw util::xml_scenario_error( msg.str() );
    }
}
double IRSAnophelesParams::RelativeAttractiveness::relativeAttractiveness(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double relAvail = exp( lPF*insecticideComponent );
    assert( relAvail>=0.0 );
    return relAvail;
}
double IRSAnophelesParams::SurvivalFactor::survivalFactor(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double killingEffect = BF + PF*insecticideComponent;
    double survivalFactor = (1.0 - killingEffect) * invBaseSurvival;
    assert( survivalFactor>=0.0 );
    return survivalFactor;
}


// ———  per-human data  ———
IRS::IRS (const TransmissionModel& tm) :
    initialInsecticide( numeric_limits<double>::signaling_NaN() )
{
    //TODO: we shouldn't really have IRS data (this class) if there's no vector
    // model, should we? Allocate dynamically or based on model?
    const VectorTransmission* vt = dynamic_cast<const VectorTransmission*>(&tm);
    if( vt != 0 ){
        const IRSParams& params = vt->getIRSParams();
        if( params.insecticideDecay.get() == 0 )
            return;     // no IRS
        // Varience factor of decay is sampled once per human: human is assumed
        // to account for most variance.
        if( params.simpleModel ){
            insecticideDecayHet = params.insecticideDecay->hetSample();
        }else{
            insecticideDecayHet = params.insecticideDecay->hetSample();
        }
    }
}

void IRS::deploy(const IRSParams& params) {
    deployTime = TimeStep::simulation;
    if( !params.simpleModel ){
        // this is sampled independently: initial insecticide content doesn't depend on handling
        initialInsecticide = params.initialInsecticide.sample();
        if( initialInsecticide < 0.0 )
            initialInsecticide = 0.0;	// avoid negative samples
        if( initialInsecticide > params.maxInsecticide )
            initialInsecticide = params.maxInsecticide;
    }
}

double IRS::relativeAttractiveness(const IRSAnophelesParams& params) const{
    double effect;
    if( params.base->simpleModel ){
        effect = (1.0 - params._relativeAttractiveness.oldDeterrency() *
            getEffectSurvival(*params.base));
    }else{
        effect = params.relativeAttractiveness( getInsecticideContent(*params.base) );
    }
    return params.byProtection( effect );
}

double IRS::preprandialSurvivalFactor(const IRSAnophelesParams& params) const{
    double effect;
    if( params.base->simpleModel ){
        effect = (1.0 - params._preprandialKillingEffect.oldEffect() *
            getEffectSurvival(*params.base));
    }else{
        effect = params.preprandialSurvivalFactor( getInsecticideContent(*params.base) );
    }
    return params.byProtection( effect );
}

double IRS::postprandialSurvivalFactor(const IRSAnophelesParams& params) const{
    double effect;
    if( params.base->simpleModel ){
        effect = (1.0 - params._postprandialKillingEffect.oldEffect() *
            getEffectSurvival(*params.base));
    }else{
        effect = params.postprandialSurvivalFactor( getInsecticideContent(*params.base) );
    }
    return params.byProtection( effect );
}

} }
