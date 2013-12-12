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

#include "interventions/Vaccine.h"
#include "Host/Human.h"
#include "util/random.h"
#include "util/errors.h"
#include "schema/interventions.h"
#include "util/StreamValidator.h"
#include "Monitoring/Surveys.h"
#include <interventions/Interventions.h>

#include <limits>
#include <cmath>

namespace OM {
namespace interventions {
using namespace OM::util;

Vaccine Vaccine::types[NumVaccineTypes];
Vaccine::Types Vaccine::reportType = NumVaccineTypes;

double Vaccine::getInitialEfficacy (size_t numPrevDoses)
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

void Vaccine::initScheduleInternal( const scnXml::ContinuousList::DeploySequence& schedule ){
    if( targetAgeTStep.size() > 0 && schedule.size() > 0 ){
        throw util::xml_scenario_error( "A vaccine effect has multiple "
                "continuous deployments. No model of how these should interact is included in OpenMalaria." );
    }
    size_t numDoses = schedule.size();
    if (numDoses > 0) {
        targetAgeTStep.resize (numDoses, TimeStep(0));
        for (size_t i = 0;i < numDoses; ++i) {
            targetAgeTStep[i] = TimeStep::fromYears( schedule[i].getTargetAgeYrs() );
        }
    }
}

void Vaccine::verifyEnabledForR_0 (){
    if( !types[PEV].active || !types[TBV].active )
        throw util::xml_scenario_error("PEV and TBV vaccines must have a "
                "description to use the insertR_0Case intervention");
}

void Vaccine::initVaccine (const scnXml::VaccineDescription& vd, Types type)
{
    if( active )
        throw util::unimplemented_exception( "multiple vaccine interventions for the same type of vaccine" );
    active = true;
    
    if( reportType == NumVaccineTypes /* the initial value */ ){
        // set to the first type described
        reportType = type;
    }

    // set efficacyB:
    efficacyB = vd.getEfficacyB().getValue();

    // set initialMeanEfficacy:
    const scnXml::VaccineDescription::InitialEfficacySequence ies = vd.getInitialEfficacy();
    initialMeanEfficacy.resize (ies.size());
    for (size_t i = 0; i < initialMeanEfficacy.size(); ++i)
        initialMeanEfficacy[i] = ies[i].getValue();

    decayFunc = DecayFunction::makeObject( vd.getDecay(), "decay" );
}

PerHumanVaccine::PerHumanVaccine(){
    for( size_t i = 0; i < Vaccine::NumVaccineTypes; ++i ){
        types[i] = PerEffectPerHumanVaccine( static_cast<Vaccine::Types>( i ) );
    }
}

PerEffectPerHumanVaccine::PerEffectPerHumanVaccine() :
    numDosesAdministered(-1),
    initialEfficacy( std::numeric_limits<double>::signaling_NaN() )
{
}

PerEffectPerHumanVaccine::PerEffectPerHumanVaccine( Vaccine::Types type ) :
        numDosesAdministered(0), initialEfficacy(0.0)
{
    if (Vaccine::types[type].active)
        hetSample = Vaccine::types[type].decayFunc->hetSample();
}

double PerEffectPerHumanVaccine::getEfficacy( Vaccine::Types type ) const{
    util::streamValidate( initialEfficacy );
    util::streamValidate( timeLastDeployment.asInt() );
    util::streamValidate( hetSample.getTMult() );
    return initialEfficacy * Vaccine::types[type].decayFunc->eval(
        TimeStep::simulation - timeLastDeployment, hetSample );
}

/// Returns true if this individual should get a vaccine dose via EPI
bool PerEffectPerHumanVaccine::getsEPIVaccination( Vaccine::Types type, TimeStep ageTSteps )const{
    const Vaccine& vacc = Vaccine::types[type];
    assert( vacc.targetAgeTStep.size() != 0 );
    // Deployment is affected by previous missed doses and mass vaccinations,
    // unlike other continuous interventions; extra test:
    return numDosesAdministered < vacc.targetAgeTStep.size()
            && vacc.targetAgeTStep[numDosesAdministered] == ageTSteps;
}

void PerEffectPerHumanVaccine::vaccinate( Deployment::Method method,
					  Vaccine::Types type ) {
    initialEfficacy = Vaccine::types[type].getInitialEfficacy(numDosesAdministered);
    util::streamValidate(initialEfficacy);

    ++numDosesAdministered;
    timeLastDeployment = TimeStep::simulation;
}

void PerHumanVaccine::vaccinate( const Host::Human& human,
                                 Deployment::Method method, Vaccine::Types type )
{
    if( method == Deployment::TIMED ||
        types[type].getsEPIVaccination( type, human.getAgeInTimeSteps() ) )
    {
        //TODO: eventually it should be unnecessary to pass type here
        types[type].vaccinate( method, type );
        
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
}
