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

#include "interventions/Vaccine.h"
#include "interventions/HumanComponents.h"
#include "Host/Human.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "schema/interventions.h"
#include "util/StreamValidator.h"
#include "Monitoring/Survey.h"

#include <limits>
#include <cmath>

namespace OM {
namespace interventions {
using namespace OM::util;
using namespace Monitoring;

vector<VaccineComponent*> VaccineComponent::params;
ComponentId VaccineComponent::reportComponent = ComponentId_pop;

VaccineComponent::VaccineComponent( ComponentId component, const scnXml::VaccineDescription& vd, Vaccine::Types type ) :
        HumanInterventionComponent(component, Report::MI_VACCINATION_CTS, Report::MI_VACCINATION_TIMED),
        type(type),
        decayFunc(DecayFunction::makeObject( vd.getDecay(), "decay" )),
        efficacyB(vd.getEfficacyB().getValue())
{
    if( type == Vaccine::BSV && ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) )
        throw util::unimplemented_exception( "blood stage vaccines (BSV) cannot be used with vivax model" );
    
    if( reportComponent == ComponentId_pop /* the initial value */ ){
        // set to the first type described
        reportComponent = component;
    }

    const scnXml::VaccineDescription::InitialEfficacySequence ies = vd.getInitialEfficacy();
    initialMeanEfficacy.resize (ies.size());
    for (size_t i = 0; i < initialMeanEfficacy.size(); ++i)
        initialMeanEfficacy[i] = ies[i].getValue();
    
    if( params.size() <= component.id ) params.resize( component.id + 1 );
    assert( params[component.id] == 0 );
    params[component.id] = this;
}

void VaccineComponent::deploy(Host::Human& human, Deployment::Method method, VaccineLimits vaccLimits) const
{
    bool administered = human.getVaccine().possiblyVaccinate( human, id(), vaccLimits );
    if( administered && VaccineComponent::reportComponent == id() ){
        Survey::current().addInt( reportMeasure(method), human, 1 );
    }
}

Component::Type VaccineComponent::componentType() const
{
    if( type == Vaccine::PEV ) return Component::PEV;
    else if( type == Vaccine::BSV ) return Component::BSV;
    else if( type == Vaccine::TBV ) return Component::TBV;
    else throw SWITCH_DEFAULT_EXCEPTION;
}

#ifdef WITHOUT_BOINC
void VaccineComponent::print_details(ostream& out) const
{
        out << id().id << "\t" << (type == Vaccine::PEV ? "PEV" : (type == Vaccine::BSV ? "BSV" : "TBV"));
}
#endif

double VaccineComponent::getInitialEfficacy (size_t numPrevDoses) const
{
    /* If initialMeanEfficacy.size or more doses have already been given, use
     * the last efficacy. */
    if (numPrevDoses >= initialMeanEfficacy.size())
        numPrevDoses = initialMeanEfficacy.size() - 1;
    double ime = initialMeanEfficacy[numPrevDoses];
    //NOTE(validation): With extra valiadation in random, the first difference is noticed here:
    util::streamValidate(ime);
    util::streamValidate(efficacyB);
    if (ime == 0.0){
        return 0.0;
    } else if (ime < 1.0) {
        double result = random::betaWithMean (ime, efficacyB);
        //NOTE(validation):: Without extra validation in random, the first difference is noticed here:
        util::streamValidate(result);
        //TODO(validation):: Why the difference? Bug in/limitation of StreamValidatior?
        // Extra memory allocation due to the extra logging causes some bad
        // memory usage to manifest differently? Am I forgetting to initialise
        // some variable?
        return result;
    } else {
        return 1.0;
    }
}

#if 0
/* R_0: this used to use _the_ vaccine configuration; this is just a check.
 * If reimplementing, use a new vaccine instance; implement some hook to ensure
 * that TBV is given to everyone after enabling the intervention, etc. */
void Vaccine::verifyEnabledForR_0 (){
    if( !types[PEV].active || !types[TBV].active )
        throw util::xml_scenario_error("PEV and TBV vaccines must have a "
                "description to use the insertR_0Case intervention");
}
#endif

// this is only used for checkpointing; loaded values are unimportant
PerEffectPerHumanVaccine::PerEffectPerHumanVaccine() :
    component(ComponentId_pop),
    numDosesAdministered(0),
    initialEfficacy( std::numeric_limits<double>::signaling_NaN() )
{
}

PerEffectPerHumanVaccine::PerEffectPerHumanVaccine( ComponentId id, const VaccineComponent& params ) :
    component( id ), numDosesAdministered(0), initialEfficacy(0.0)
{
    hetSample = params.decayFunc->hetSample();
}

double PerHumanVaccine::getFactor( Vaccine::Types type ) const{
    double factor = 1.0;
    for( EffectList::const_iterator effect = effects.begin(); effect != effects.end(); ++effect ){
        if( VaccineComponent::getParams(effect->component).type == type ){
            SimTime age = sim::ts1() - effect->timeLastDeployment;  // implies age 1 TS on first use
            double decayFactor = VaccineComponent::getParams(effect->component)
                .decayFunc->eval( age, effect->hetSample );
            factor *= 1.0 - effect->initialEfficacy * decayFactor;
        }
    }
    return factor;
}

bool PerHumanVaccine::possiblyVaccinate( const Host::Human& human,
        ComponentId componentId, interventions::VaccineLimits vaccLimits )
{
    PerEffectPerHumanVaccine* effect = 0;
    for( EffectList::iterator it = effects.begin(); it != effects.end(); ++it ){
        if( it->component == componentId ){
            effect = &*it;
            break;
        }
    }
    
    uint32_t numDosesAdministered = (effect == 0) ? 0 : effect->numDosesAdministered;
    if( numDosesAdministered < vaccLimits.minPrevDoses ||
        numDosesAdministered >= vaccLimits.maxCumDoses )
        return false;   // no vaccination (this replaces the old schedule for continuous doses)
    
    const VaccineComponent& params = VaccineComponent::getParams( componentId );
    
    if( effect == 0 ){
        effects.push_back( PerEffectPerHumanVaccine( componentId, params ) );
        effect = &effects.back();
    }
    
    effect->initialEfficacy = params.getInitialEfficacy(numDosesAdministered);
    util::streamValidate(effect->initialEfficacy);
    
    effect->numDosesAdministered = numDosesAdministered + 1;
    effect->timeLastDeployment = sim::nowOrTs1();
    
    return true;
}

}
}
