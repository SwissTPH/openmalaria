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

Vaccine* Vaccine::params[NumVaccineTypes];
Vaccine::Types Vaccine::reportType = NumVaccineTypes;

Vaccine::Vaccine(const scnXml::VaccineDescription& vd, Vaccine::Types type ) :
        decayFunc(DecayFunction::makeObject( vd.getDecay(), "decay" )),
        efficacyB(vd.getEfficacyB().getValue())
{
    if( params[type] != 0 )
        throw util::unimplemented_exception( "multiple vaccine interventions for the same type of vaccine" );
    if( type == BSV && ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) )
        throw util::unimplemented_exception( "blood stage vaccines (BSV) cannot be used with vivax model" );
    
    if( reportType == NumVaccineTypes /* the initial value */ ){
        // set to the first type described
        reportType = type;
    }

    const scnXml::VaccineDescription::InitialEfficacySequence ies = vd.getInitialEfficacy();
    initialMeanEfficacy.resize (ies.size());
    for (size_t i = 0; i < initialMeanEfficacy.size(); ++i)
        initialMeanEfficacy[i] = ies[i].getValue();
    
    params[type] = this;
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

PerEffectPerHumanVaccine::PerEffectPerHumanVaccine() :
    numDosesAdministered(0),
    initialEfficacy( std::numeric_limits<double>::signaling_NaN() )
{
}

PerEffectPerHumanVaccine::PerEffectPerHumanVaccine( const Vaccine& params ) :
        numDosesAdministered(0), initialEfficacy(0.0)
{
    hetSample = params.decayFunc->hetSample();
}

double PerHumanVaccine::getFactor( Vaccine::Types type ) const{
    if( types[type] == 0 ) return 1.0;  // never deployed
    const PerEffectPerHumanVaccine& effect = *types[type];
    TimeStep age = TimeStep::simulation - effect.timeLastDeployment;
    double decayFactor = Vaccine::getParams(type).decayFunc->eval( age, effect.hetSample );
    return 1.0 - effect.initialEfficacy * decayFactor;
}

void PerHumanVaccine::possiblyVaccinate( const Host::Human& human,
                                 Deployment::Method method, Vaccine::Types type,
                                 interventions::VaccineLimits vaccLimits )
{
    PerEffectPerHumanVaccine* effect = types[type];
    
    uint32_t numDosesAdministered = (effect == 0) ? 0 : effect->numDosesAdministered;
    if( numDosesAdministered < vaccLimits.minPrevDoses ||
        numDosesAdministered >= vaccLimits.maxCumDoses )
        return;         // no vaccination (this replaces the old schedule for continuous doses)
    
    const Vaccine& params = Vaccine::getParams( type );
    
    if( effect == 0 ){
        types[type] = new PerEffectPerHumanVaccine( params );
        effect = types[type];
    }
    
    effect->initialEfficacy = params.getInitialEfficacy(numDosesAdministered);
    util::streamValidate(effect->initialEfficacy);
    
    effect->numDosesAdministered = numDosesAdministered + 1;
    effect->timeLastDeployment = TimeStep::simulation;
    
    if( Vaccine::reportType == type ){
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
