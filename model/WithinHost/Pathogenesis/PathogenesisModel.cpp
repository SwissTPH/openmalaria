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

#include "WithinHost/Pathogenesis/Submodels.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include <Parameters.h>
#include <schema/healthSystem.h>

#include <cmath>
using namespace std;

namespace OM { namespace WithinHost { namespace Pathogenesis {

using namespace OM::util;

//BEGIN static
double PathogenesisModel::indirRiskCoFactor_18;
double PathogenesisModel::sevMal_21;
double PathogenesisModel::comorbintercept_24;
double PathogenesisModel::critAgeComorb_30;

AgeGroupInterpolation* PathogenesisModel::NMF_incidence = AgeGroupInterpolation::dummyObject();

bool opt_predetermined_episodes = false, opt_mueller_pres_model = false;


void PathogenesisModel::init( const Parameters& parameters, const scnXml::Clinical& clinical ) {
    indirRiskCoFactor_18=(1-exp(-parameters[Parameters::INDIRECT_RISK_COFACTOR]));
    sevMal_21=parameters[Parameters::SEVERE_MALARIA_THRESHHOLD];
    comorbintercept_24=1-exp(-parameters[Parameters::COMORBIDITY_INTERCEPT]);
    critAgeComorb_30=parameters[Parameters::CRITICAL_AGE_FOR_COMORBIDITY];

    if (util::ModelOptions::option (util::PREDETERMINED_EPISODES)) {
        opt_predetermined_episodes = true;
        //no separate init:
        PyrogenPathogenesis::init( parameters );
    } else {
        if (util::ModelOptions::option (util::MUELLER_PRESENTATION_MODEL)){
            opt_mueller_pres_model = true;
            MuellerPathogenesis::init( parameters );
        }else{
            PyrogenPathogenesis::init( parameters );
        }
    }
    
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !clinical.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc = clinical.getNonMalariaFevers().get();
        NMF_incidence = AgeGroupInterpolation::makeObject( nmfDesc.getIncidence(), "incidence" );
    }
}
void PathogenesisModel::cleanup() {
    AgeGroupInterpolation::freeObject( NMF_incidence );
}

PathogenesisModel* PathogenesisModel::createPathogenesisModel(double cF) {
    if (opt_predetermined_episodes) {
        return new PredetPathogenesis(cF);
    }
    else {
        if (opt_mueller_pres_model) {
            return new MuellerPathogenesis(cF);
        }
        else {
            return new PyrogenPathogenesis(cF);
        }
    }
}
//END static


PathogenesisModel::PathogenesisModel(double cF) :
        _comorbidityFactor(cF)
{}

Pathogenesis::StatePair PathogenesisModel::determineState (double ageYears, double timeStepMaxDensity, double endDensity) {
    double prEpisode = getPEpisode(timeStepMaxDensity, endDensity);
    StatePair result;

    //Decide whether a clinical episode occurs and if so, which type
    if ((random::uniform_01()) < prEpisode) {
        //Fixed severe threshold
        double severeMalThreshold=sevMal_21+1;
        double prSevereEpisode=1-1/(1+timeStepMaxDensity/severeMalThreshold);

        if (random::uniform_01() < prSevereEpisode)
            result.state = STATE_SEVERE;
        else {
            double pCoinfection=comorbintercept_24/(1+ageYears/critAgeComorb_30);
            pCoinfection*=_comorbidityFactor;

            if (random::uniform_01() < pCoinfection)
                result.state = STATE_COINFECTION;
            else
                result.state = STATE_MALARIA;
        }

        /* Indirect mortality
           IndirectRisk is the probability of dying from indirect effects of malaria
           conditional on not having an acute attack of malaria
        */
        double indirectRisk=indirRiskCoFactor_18/(1+ageYears/critAgeComorb_30);
        indirectRisk*=_comorbidityFactor;
        if (random::uniform_01() < indirectRisk){
            result.indirectMortality = true;
        }
    } else if ( NMF_incidence->isSet() ) {
        if( random::uniform_01() < NMF_incidence->eval( ageYears ) ){
            result.state = STATE_NMF;
        }
    }
    return result;
}


void PathogenesisModel::checkpoint (istream& stream) {
    _comorbidityFactor & stream;
}
void PathogenesisModel::checkpoint (ostream& stream) {
    _comorbidityFactor & stream;
}

} } }
