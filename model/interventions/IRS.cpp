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

#include "interventions/IRS.h"
#include "Host/Human.h"
#include "util/random.h"
#include "util/errors.h"
#include "Monitoring/Surveys.h"
#include "util/SpeciesIndexChecker.h"

#include "R_nmath/qnorm.h"
#include <cmath>

namespace OM { namespace interventions {
    using util::random::poisson;

vector<IRSEffect*> IRSEffect::effectsByIndex;

IRSEffect::IRSEffect( size_t index, const scnXml::IRSDescription& elt,
        const map<string,size_t>& species_name_map ) :
        Transmission::HumanVectorInterventionEffect(index)
{
    initialInsecticide.setParams( elt.getInitialInsecticide() );
    const double maxProp = 0.999;       //NOTE: this could be exposed in XML, but probably doesn't need to be
    maxInsecticide = R::qnorm5(maxProp, initialInsecticide.getMu(), initialInsecticide.getSigma(), true, false);
    insecticideDecay = DecayFunction::makeObject( elt.getInsecticideDecay(),
                                                  "interventions.human.IRS.description.insecticideDecay" );
    
    typedef scnXml::IRSDescription::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    species.resize( species_name_map.size() );
    util::SpeciesIndexChecker checker( "IRS intervention", species_name_map );
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ) {
        species[checker.getIndex(it->getMosquito())].init (*it, maxInsecticide);
    }
    checker.checkNoneMissed();
    
    if( effectsByIndex.size() <= index ) effectsByIndex.resize( index+1, 0 );
    effectsByIndex[index] = this;
}

void IRSEffect::deploy( Host::Human& human, Deployment::Method method )const{
    human.perHostTransmission.deployEffect(*this);
    if( method == interventions::Deployment::TIMED ){
        Monitoring::Surveys.getSurvey(human.isInAnyCohort()).reportMassIRS( human.getMonitoringAgeGroup(), 1 );
    }else if( method == interventions::Deployment::CTS ){
        //TODO(monitoring): report
    }else throw SWITCH_DEFAULT_EXCEPTION;
}

Effect::Type IRSEffect::effectType()const{ return Effect::IRS; }

PerHostInterventionData* IRSEffect::makeHumanPart() const{
    return new HumanIRS( *this );
}
PerHostInterventionData* IRSEffect::makeHumanPart( istream& stream, size_t index ) const{
    return new HumanIRS( stream, index );
}

void IRSEffect::IRSAnopheles::init(const scnXml::IRSDescription::AnophelesParamsType& elt,
                                   double maxInsecticide)
{
    _relativeAttractiveness.init( elt.getDeterrency() );
    _preprandialKillingEffect.init( elt.getPreprandialKillingEffect(), false, maxInsecticide );
    _postprandialKillingEffect.init( elt.getPostprandialKillingEffect(), true, maxInsecticide );
    // Simpler version of ITN usage/action:
    double propActive = elt.getPropActive();
    assert( propActive >= 0.0 && propActive <= 1.0 );
    proportionProtected = propActive;
    proportionUnprotected = 1.0 - proportionProtected;
}

inline bool inRange01( double x ){
    return x>=0.0 && x<= 1.0;
}
IRSEffect::IRSAnopheles::RelativeAttractiveness::RelativeAttractiveness() :
    lPF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() )
{}
void IRSEffect::IRSAnopheles::RelativeAttractiveness::init(const scnXml::IRSDeterrency& elt){
    double PF = elt.getInsecticideFactor();
    insecticideScaling = elt.getInsecticideScalingFactor();
    if( !( PF > 0.0) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams.relativeAvailability: expected insecticideFactor to be positive.";
        //TODO: These constraints were required. But they're too strong.
        // Now need to work out which should still be imposed.
        cerr << msg.str() << endl;
        //throw util::xml_scenario_error( msg.str() );
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
    if( lPF == lPF ){
        throw util::unimplemented_exception( "multiple IRS interventions" );
    }
    lPF = log( PF );
}
IRSEffect::IRSAnopheles::SurvivalFactor::SurvivalFactor() :
    BF( numeric_limits< double >::signaling_NaN() ),
    PF( numeric_limits< double >::signaling_NaN() ),
    insecticideScaling( numeric_limits< double >::signaling_NaN() ),
    invBaseSurvival( numeric_limits< double >::signaling_NaN() )
{}
void IRSEffect::IRSAnopheles::SurvivalFactor::init(const scnXml::IRSKillingEffect& elt,
                                                   bool postPrandial, double maxInsecticide){
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
    double pmax = 1.0-exp(-maxInsecticide*insecticideScaling);
    if( !( PF >= 0.0 && BF+PF*pmax <= 1.0 ) ){
        ostringstream msg;
        msg << "IRS.description.anophelesParams." << (postPrandial?"post":"pre")
            << "killingFactor: expected insecticideFactor≥0, baseFactor+"
            <<pmax<<"×insecticideFactor≤1";
        //TODO: These constraints were required. But they're too strong.
        // Now need to work out which should still be imposed.
        cerr << msg.str() << endl;
        //throw util::xml_scenario_error( msg.str() );
    }
}
double IRSEffect::IRSAnopheles::RelativeAttractiveness::relativeAttractiveness(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double relAvail = exp( lPF*insecticideComponent );
    //TODO: limits
    //assert( relAvail>=0.0 );
    if (relAvail < 0.0)
        return 0.0;
    return relAvail;
}
double IRSEffect::IRSAnopheles::SurvivalFactor::survivalFactor(
        double insecticideContent
) const {
    double insecticideComponent = 1.0 - exp(-insecticideContent*insecticideScaling);
    double killingEffect = BF + PF*insecticideComponent;
    double survivalFactor = (1.0 - killingEffect) * invBaseSurvival;
    //TODO: limits
    //assert( survivalFactor >= 0.0 );
    if (survivalFactor < 0.0)
        return 0.0;
    else if (survivalFactor > 1.0)
        return 1.0;
    return survivalFactor;
}

double IRSEffect::sampleInitialInsecticide() const{
    double sample = initialInsecticide.sample();
    if( sample < 0.0 )
        sample = 0.0;       // avoid negative samples
    if( sample > maxInsecticide )
        sample = maxInsecticide;
    return sample;
}


// ———  per-human data  ———
HumanIRS::HumanIRS( const IRSEffect& params ) :
    PerHostInterventionData( params.getIndex() ),
    initialInsecticide( params.sampleInitialInsecticide() )
{
    // Varience factor of decay is sampled once per human: human is assumed
    // to account for most variance.
    insecticideDecayHet = params.insecticideDecay->hetSample();
}

void HumanIRS::redeploy( const OM::Transmission::HumanVectorInterventionEffect& params ) {
    deployTime = TimeStep::simulation;
    initialInsecticide = dynamic_cast<const IRSEffect*>(&params)->sampleInitialInsecticide();
}

void HumanIRS::update(){
}

double HumanIRS::relativeAttractiveness(size_t speciesIndex) const{
    const IRSEffect& params = *IRSEffect::effectsByIndex[index];
    const IRSEffect::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.relativeAttractiveness( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}

double HumanIRS::preprandialSurvivalFactor(size_t speciesIndex) const{
    const IRSEffect& params = *IRSEffect::effectsByIndex[index];
    const IRSEffect::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.preprandialSurvivalFactor( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}

double HumanIRS::postprandialSurvivalFactor(size_t speciesIndex) const{
    const IRSEffect& params = *IRSEffect::effectsByIndex[index];
    const IRSEffect::IRSAnopheles& anoph = params.species[speciesIndex];
    double effect = anoph.postprandialSurvivalFactor( getInsecticideContent(params) );
    return anoph.byProtection( effect );
}

void HumanIRS::checkpoint( ostream& stream ){
    deployTime & stream;
    initialInsecticide & stream;
    insecticideDecayHet & stream;
}
HumanIRS::HumanIRS( istream& stream, size_t index ) : PerHostInterventionData( index )
{
    deployTime & stream;
    initialInsecticide & stream;
    insecticideDecayHet & stream;
}

} }
