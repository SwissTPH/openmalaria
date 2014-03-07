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

#include "WithinHost/WHFalciparum.h"
#include "WithinHost/DescriptiveWithinHost.h"
#include "WithinHost/CommonWithinHost.h"
#include "WithinHost/Infection/DummyInfection.h"
#include "WithinHost/Infection/EmpiricalInfection.h"
#include "WithinHost/Infection/MolineauxInfection.h"
#include "WithinHost/Infection/PennyInfection.h"
#include "WithinHost/Pathogenesis/PathogenesisModel.h"
#include "WithinHost/Diagnostic.h"
#include "WithinHost/Treatments.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "schema/scenario.h"

#include <cmath>
#include <boost/format.hpp>
#include <gsl/gsl_cdf.h>


namespace OM {
namespace WithinHost {

using namespace OM::util;

double WHFalciparum::sigma_i;
double WHFalciparum::immPenalty_22;
double WHFalciparum::asexImmRemain;
double WHFalciparum::immEffectorRemain;
int WHFalciparum::_ylagLen = 0;

// -----  static functions  -----

void WHFalciparum::init( const OM::Parameters& parameters, const scnXml::Scenario& scenario ) {
    sigma_i=sqrt(parameters[Parameters::SIGMA_I_SQ]);
    immPenalty_22=1-exp(parameters[Parameters::IMMUNITY_PENALTY]);
    immEffectorRemain=exp(-parameters[Parameters::IMMUNE_EFFECTOR_DECAY]);
    asexImmRemain=exp(-parameters[Parameters::ASEXUAL_IMMUNITY_DECAY]);
    
    _ylagLen = TimeStep::intervalsPer5Days.asInt() * 4;
    
    //NOTE: should also call cleanup() on the PathogenesisModel, but it only frees memory which the OS does anyway
    Pathogenesis::PathogenesisModel::init( parameters, scenario.getModel().getClinical() );
    
    /*
    The detection limit (in parasites/ul) is currently the same for PCR and for microscopy
    TODO: in fact the detection limit in Garki should be the same as the PCR detection limit
    The density bias allows the detection limit for microscopy to be higher for other sites
    */
    double densitybias;
    if (util::ModelOptions::option (util::GARKI_DENSITY_BIAS)) {
        densitybias = parameters[Parameters::DENSITY_BIAS_GARKI];
    } else {
        if( scenario.getAnalysisNo().present() ){
            int analysisNo = scenario.getAnalysisNo().get();
            if ((analysisNo >= 22) && (analysisNo <= 30)) {
                cerr << "Warning: these analysis numbers used to mean use Garki density bias. If you do want to use this, specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong." << endl;
            }
        }
        densitybias = parameters[Parameters::DENSITY_BIAS_NON_GARKI];
    }
    double detectionLimit=scenario.getMonitoring().getSurveys().getDetectionLimit()*densitybias;
    Diagnostic::default_.setDeterministic( detectionLimit );
    
    Infection::init( parameters, scenario.getModel().getParameters().getLatentp() );
}


// -----  Non-static  -----

WHFalciparum::WHFalciparum( double comorbidityFactor ):
    WHInterface(),
    _cumulativeh(0.0), _cumulativeY(0.0), _cumulativeYlag(0.0),
    totalDensity(0.0), timeStepMaxDensity(0.0),
    pathogenesisModel( Pathogenesis::PathogenesisModel::createPathogenesisModel( comorbidityFactor ) )
{
    _innateImmSurvFact = exp(-random::gauss(0, sigma_i));
    
    _ylag.assign (_ylagLen, 0.0);
}

WHFalciparum::~WHFalciparum()
{
}

double WHFalciparum::probTransmissionToMosquito( TimeStep ageTimeSteps, double tbvFactor ) const{
    /* This model (often referred to as the gametocyte model) was designed for
    5-day timesteps. We use the same model (sampling 10, 15 and 20 days ago)
    for 1-day timesteps to avoid having to design and analyse a new model.
    Description: AJTMH pp.32-33 */
    
    /* Note: we don't allow for treatment which clears gametocytes (e.g.
     * Primaquine). Apparently Primaquine is not commonly used in P falciparum
     * treatment, but for vivax the effect may be important. */
    
    if (ageTimeSteps.inDays() <= 20 || TimeStep::simulation.inDays() <= 20){
        // We need at least 20 days history (_ylag) to calculate infectiousness;
        // assume no infectiousness if we don't have this history.
        // Note: human not updated on DOB so age must be >20 days.
        return 0.0;
    }
    
    //Infectiousness parameters: see AJTMH p.33, tau=1/sigmag**2 
    static const double beta1=1.0;
    static const double beta2=0.46;
    static const double beta3=0.17;
    static const double tau= 0.066;
    static const double mu= -8.1;
    
    // Take weighted sum of total asexual blood stage density 10, 15 and 20 days
    // before. We have 20 days history, so use mod_nn:
    int firstIndex = TimeStep::simulation.asInt()-2*TimeStep::intervalsPer5Days.asInt() + 1;
    double x = beta1 * _ylag[mod_nn(firstIndex, _ylagLen)]
            + beta2 * _ylag[mod_nn(firstIndex-TimeStep::intervalsPer5Days.asInt(), _ylagLen)]
            + beta3 * _ylag[mod_nn(firstIndex-2*TimeStep::intervalsPer5Days.asInt(), _ylagLen)];
    if (x < 0.001){
        return 0.0;
    }
    
    double zval=(log(x)+mu)/sqrt(1.0/tau);
    double pone = gsl_cdf_ugaussian_P(zval);
    double transmit=(pone*pone);
    //transmit has to be between 0 and 1
    transmit=std::max(transmit, 0.0);
    transmit=std::min(transmit, 1.0);
    
    //    Include here the effect of transmission-blocking vaccination
    double result = transmit * tbvFactor;
    util::streamValidate( result );
    return result;
}

bool WHFalciparum::diagnosticDefault() const{
    return Diagnostic::default_.isPositive( totalDensity );
}

void WHFalciparum::treatment(TreatmentId treatId){
    const Treatments& treat = Treatments::select( treatId );
    for( vector<Treatments::Action>::const_iterator it =
        treat.getEffects().begin(), end = treat.getEffects().end();
        it != end; ++it )
    {
        if( it->timesteps == TimeStep(-1) ){
            // act immediately
            //TODO: which timestep to measure "is blood stage" by?
            clearInfections( it->stage );
        }else{
            switch( it->stage ){
                case Treatments::BOTH:
                    treatExpiryLiver = max( treatExpiryLiver, TimeStep::simulation + it->timesteps );
                    // don't break; also do blood below:
                case Treatments::BLOOD:
                    treatExpiryBlood = max( treatExpiryBlood, TimeStep::simulation + it->timesteps );
                    break;
                case Treatments::LIVER:
                    treatExpiryLiver = max( treatExpiryLiver, TimeStep::simulation + it->timesteps );
                    break;
                case Treatments::NONE:
                    /*do nothing*/;
            }
        }
    }
}

Pathogenesis::StatePair WHFalciparum::determineMorbidity(double ageYears){
    Pathogenesis::StatePair result =
            pathogenesisModel->determineState( ageYears, timeStepMaxDensity, totalDensity );
    
    /* Note: this model can easily be re-enabled, but is not used and not considered to be a good model.
    if( (result.state & Pathogenesis::MALARIA) && util::ModelOptions::option( util::PENALISATION_EPISODES ) ){
        // This does immunity penalisation:
        _cumulativeY = _cumulativeYlag - immPenalty_22*(_cumulativeY-_cumulativeYlag);
        if (_cumulativeY < 0) {
            _cumulativeY=0.0;
        }
    }*/
    
    return result;
}


// -----  immunity  -----

void WHFalciparum::updateImmuneStatus() {
    if (immEffectorRemain < 1) {
        _cumulativeh*=immEffectorRemain;
        _cumulativeY*=immEffectorRemain;
    }
    if (asexImmRemain < 1) {
        _cumulativeh*=asexImmRemain/
                      (1+(_cumulativeh*(1-asexImmRemain) * Infection::invCumulativeHstar));
        _cumulativeY*=asexImmRemain/
                      (1+(_cumulativeY*(1-asexImmRemain) * Infection::invCumulativeYstar));
    }
    _cumulativeYlag = _cumulativeY;
}


// -----  Summarize  -----

bool WHFalciparum::summarize (Monitoring::Survey& survey, Monitoring::AgeGroup ageGroup) {
    pathogenesisModel->summarize( survey, ageGroup );
    InfectionCount count = countInfections();
    if (count.total != 0) {
        survey.reportInfectedHosts(ageGroup,1);
        survey.addToInfections(ageGroup, count.total);
        survey.addToPatentInfections(ageGroup, count.patent);
    }
    // Treatments in the old ImmediateOutcomes clinical model clear infections immediately
    // (and are applied after update()); here we report the last calculated density.
    if (diagnosticDefault()) {
        survey.reportPatentHosts(ageGroup, 1);
        survey.addToLogDensity(ageGroup, log(totalDensity));
        return true;
    }
    return false;
}


void WHFalciparum::checkpoint (istream& stream) {
    WHInterface::checkpoint( stream );
    _innateImmSurvFact & stream;
    _cumulativeh & stream;
    _cumulativeY & stream;
    _cumulativeYlag & stream;
    totalDensity & stream;
    timeStepMaxDensity & stream;
    _ylag & stream;
    (*pathogenesisModel) & stream;
}
void WHFalciparum::checkpoint (ostream& stream) {
    WHInterface::checkpoint( stream );
    _innateImmSurvFact & stream;
    _cumulativeh & stream;
    _cumulativeY & stream;
    _cumulativeYlag & stream;
    totalDensity & stream;
    timeStepMaxDensity & stream;
    _ylag & stream;
    (*pathogenesisModel) & stream;
}

}
}
