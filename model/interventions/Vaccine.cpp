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
#include "Host/Human.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "schema/interventions.h"
#include "util/StreamValidator.h"
#include "Monitoring/Surveys.h"

#include <limits>
#include <cmath>

namespace OM {
namespace interventions {
using namespace OM::util;

vector<Vaccine*> Vaccine::params;
EffectId Vaccine::reportEffect = EffectId_pop;

Vaccine::Vaccine(const scnXml::VaccineDescription& vd, Vaccine::Types type, EffectId effect) :
        type(type),
        decayFunc(DecayFunction::makeObject( vd.getDecay(), "decay" )),
        efficacyB(vd.getEfficacyB().getValue())
{
    if( type == BSV && ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) )
        throw util::unimplemented_exception( "blood stage vaccines (BSV) cannot be used with vivax model" );
    
    if( reportEffect == EffectId_pop /* the initial value */ ){
        // set to the first type described
        reportEffect = effect;
    }

    const scnXml::VaccineDescription::InitialEfficacySequence ies = vd.getInitialEfficacy();
    initialMeanEfficacy.resize (ies.size());
    for (size_t i = 0; i < initialMeanEfficacy.size(); ++i)
        initialMeanEfficacy[i] = ies[i].getValue();
    
    if( params.size() <= effect.id ) params.resize( effect.id + 1 );
    assert( params[effect.id] == 0 );
    params[effect.id] = this;
}

double Vaccine::getInitialEfficacy (size_t numPrevDoses) const
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
    effect(EffectId_pop),
    numDosesAdministered(0),
    initialEfficacy( std::numeric_limits<double>::signaling_NaN() )
{
}

PerEffectPerHumanVaccine::PerEffectPerHumanVaccine( EffectId id, const Vaccine& params ) :
    effect( id ), numDosesAdministered(0), initialEfficacy(0.0)
{
    hetSample = params.decayFunc->hetSample();
}

double PerHumanVaccine::getFactor( Vaccine::Types type ) const{
    double factor = 1.0;
    for( EffectList::const_iterator effect = effects.begin(); effect != effects.end(); ++effect ){
        if( Vaccine::getParams(effect->effect).type == type ){
            TimeStep age = TimeStep::simulation - effect->timeLastDeployment;
            double decayFactor = Vaccine::getParams(effect->effect).decayFunc->eval( age, effect->hetSample );
            factor *= 1.0 - effect->initialEfficacy * decayFactor;
        }
    }
    return factor;
}

void PerHumanVaccine::possiblyVaccinate( const Host::Human& human,
                                 Deployment::Method method, EffectId effectId,
                                 interventions::VaccineLimits vaccLimits )
{
    PerEffectPerHumanVaccine* effect = 0;
    for( EffectList::iterator it = effects.begin(); it != effects.end(); ++it ){
        if( it->effect == effectId ){
            effect = &*it;
            break;
        }
    }
    
    uint32_t numDosesAdministered = (effect == 0) ? 0 : effect->numDosesAdministered;
    if( numDosesAdministered < vaccLimits.minPrevDoses ||
        numDosesAdministered >= vaccLimits.maxCumDoses )
        return;         // no vaccination (this replaces the old schedule for continuous doses)
    
    const Vaccine& params = Vaccine::getParams( effectId );
    
    if( effect == 0 ){
        effects.push_back( PerEffectPerHumanVaccine( effectId, params ) );
        effect = &effects.back();
    }
    
    effect->initialEfficacy = params.getInitialEfficacy(numDosesAdministered);
    util::streamValidate(effect->initialEfficacy);
    
    effect->numDosesAdministered = numDosesAdministered + 1;
    effect->timeLastDeployment = TimeStep::simulation;
    
    if( Vaccine::reportEffect == effectId ){
        if( method == Deployment::TIMED )
            Monitoring::Surveys.getSurvey(human.isInAnyCohort())
                .reportMassVaccinations (human.getMonitoringAgeGroup(), 1);
        else
            Monitoring::Surveys.getSurvey(human.isInAnyCohort())
                .reportEPIVaccinations (human.getMonitoringAgeGroup(), 1);
    }
}

}
}
