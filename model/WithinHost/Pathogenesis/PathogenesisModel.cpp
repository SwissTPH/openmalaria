/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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
#include "util/AgeGroupInterpolation.h"
#include "Parameters.h"
#include "schema/healthSystem.h"

#include <cmath>
using namespace std;

namespace OM { namespace WithinHost { namespace Pathogenesis {

using namespace OM::util;
using util::AgeGroupInterpolator;

//BEGIN static
/// Comorbidity prevalence at birth as a risk factor for indirect mortality
double pg_indirRiskCoFactor;
/// sevMal: critical density for severe malaria bout (Y*B1)
double pg_severeMalThreshold;
/// Comorbidity prevalence at birth as a risk factor for severe
double pg_comorbIntercept;
/// One over critical age for co-morbidity (for both severe and indirect)
double pg_inv_critAgeComorb;

/// Rate of Non-Malaria Fever incidence by age. Non-seasonal.
AgeGroupInterpolator pg_NMF_incidence;

bool opt_predetermined_episodes = false, opt_mueller_pres_model = false;


void PathogenesisModel::init( const Parameters& parameters, const scnXml::Clinical& clinical, bool nmfOnly ){
    if( util::ModelOptions::option( util::NON_MALARIA_FEVERS ) ){
        if( !clinical.getNonMalariaFevers().present() ){
            throw util::xml_scenario_error("NonMalariaFevers element of model->clinical required");
        }
        const scnXml::Clinical::NonMalariaFeversType& nmfDesc = clinical.getNonMalariaFevers().get();
        pg_NMF_incidence.set( nmfDesc.getIncidence(), "incidence" );
    }
    if( nmfOnly ) return;
    
    pg_indirRiskCoFactor = 1 - exp(-parameters[Parameters::INDIRECT_RISK_COFACTOR]);
    pg_severeMalThreshold = parameters[Parameters::SEVERE_MALARIA_THRESHHOLD] + 1;
    pg_comorbIntercept = 1 - exp(-parameters[Parameters::COMORBIDITY_INTERCEPT]);
    pg_inv_critAgeComorb = 1 / parameters[Parameters::CRITICAL_AGE_FOR_COMORBIDITY];

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
    double pMalariaFever = getPEpisode(timeStepMaxDensity, endDensity);
    StatePair result;
    //TODO(performance): would using a single RNG sample and manipulating probabilities be faster?
    //Decide whether a clinical episode occurs and if so, which type
    if( random::bernoulli( pMalariaFever ) ){
        double prSevereEpisode = timeStepMaxDensity / (timeStepMaxDensity + pg_severeMalThreshold);
        double comorb_factor = _comorbidityFactor / (1.0 + ageYears * pg_inv_critAgeComorb);
        
        if( random::bernoulli( prSevereEpisode ) )
            result.state = STATE_SEVERE;
        else {
            double pCoinfection = pg_comorbIntercept * comorb_factor;
            if( random::bernoulli( pCoinfection ) )
                result.state = STATE_COINFECTION;
            else
                result.state = STATE_MALARIA;
        }

        // Indirect mortality:
        // IndirectRisk is the probability of dying from the indirect effects
        // of malaria conditional on not having an acute attack of malaria
        double indirectRisk = pg_indirRiskCoFactor * comorb_factor;
        if( random::bernoulli( indirectRisk ) )
            result.indirectMortality = true;        
    }else{
        result.state = sampleNMF( ageYears );
    }
    return result;
}

Pathogenesis::State PathogenesisModel::sampleNMF( double ageYears ){
    if ( pg_NMF_incidence.isSet() ) {
        double pNMF = pg_NMF_incidence.eval( ageYears );
        if( random::bernoulli( pNMF ) )
            return Pathogenesis::STATE_NMF;
    }
    return Pathogenesis::NONE;
}


void PathogenesisModel::checkpoint (istream& stream) {
    _comorbidityFactor & stream;
}
void PathogenesisModel::checkpoint (ostream& stream) {
    _comorbidityFactor & stream;
}

} } }
