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

#include "interventions/GVI.h"
#include "Host/Human.h"
#include "util/SpeciesIndexChecker.h"
#include <cmath>

namespace OM { namespace interventions {

GVIParams::GVIParams( size_t index, const scnXml::GVIDescription& elt,
        const map<string,size_t>& species_name_map ) : HumanVectorInterventionParams(index)
{
    decay = DecayFunction::makeObject( elt.getDecay(), "interventions.human.vector.decay" );
    
    typedef scnXml::GVIDescription::AnophelesParamsSequence AP;
    const AP& ap = elt.getAnophelesParams();
    species.resize( species_name_map.size() );
    util::SpeciesIndexChecker checker( "GVI intervention", species_name_map );
    for( AP::const_iterator it = ap.begin(); it != ap.end(); ++it ) {
        species[checker.getIndex(it->getMosquito())].init (*it);
    }
    checker.checkNoneMissed();
}

 void GVIParams::deploy( Host::Human& human, Deployment::Method method )const{
     human.perHostTransmission.interventions.deploy(*this);
     human.reportDeployment( Effect::GVI, method );
 }
 
 HumanVectorIntervention* GVIParams::makeHumanPart() const{
     return new HumanGVI( *this );
 }
 HumanVectorIntervention* GVIParams::makeHumanPart( istream& stream, size_t index ) const{
     return new HumanGVI( stream, index );
 }

void GVIParams::GVIAnopheles::init(const scnXml::GVIDescription::AnophelesParamsType& elt)
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
HumanGVI::HumanGVI ( const GVIParams& params ) :
    HumanVectorIntervention( params.getIndex() ),
    initialInsecticide( 0.0 )   // start with no insecticide (for monitoring)
{
    // Varience factor of decay is sampled once per human: human is assumed
    // to account for most variance.
    decayHet = params.decay->hetSample();
}

void HumanGVI::deploy( const HumanVectorInterventionParams& params ) {
    deployTime = TimeStep::simulation;
}

double HumanGVI::relativeAttractiveness(const HumanInterventionEffect& gen_params, size_t speciesIndex) const{
    assert( dynamic_cast<const GVIParams*>(&gen_params) != 0 );
    const GVIParams& params = *dynamic_cast<const GVIParams*>(&gen_params);
    const GVIParams::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph._relativeAttractiveness *
            getEffectSurvival(params));
    return anoph.byProtection( effect );
}

double HumanGVI::preprandialSurvivalFactor(const HumanInterventionEffect& gen_params, size_t speciesIndex) const{
    assert( dynamic_cast<const GVIParams*>(&gen_params) != 0 );
    const GVIParams& params = *dynamic_cast<const GVIParams*>(&gen_params);
    const GVIParams::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph._preprandialKillingEffect *
            getEffectSurvival(params));
    return anoph.byProtection( effect );
}

double HumanGVI::postprandialSurvivalFactor(const HumanInterventionEffect& gen_params, size_t speciesIndex) const{
    assert( dynamic_cast<const GVIParams*>(&gen_params) != 0 );
    const GVIParams& params = *dynamic_cast<const GVIParams*>(&gen_params);
    const GVIParams::GVIAnopheles& anoph = params.species[speciesIndex];
    double effect = (1.0 - anoph._postprandialKillingEffect *
            getEffectSurvival(params));
    return anoph.byProtection( effect );
}

void HumanGVI::checkpoint( ostream& stream ){
    deployTime & stream;
    initialInsecticide & stream;
    decayHet & stream;
}
HumanGVI::HumanGVI( istream& stream, size_t index ) : HumanVectorIntervention( index )
{
    deployTime & stream;
    initialInsecticide & stream;
    decayHet & stream;
}

} }
