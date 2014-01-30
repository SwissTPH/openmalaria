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
#include "WithinHost/DescriptiveIPTWithinHost.h"
#include "WithinHost/CommonWithinHost.h"
#include "WithinHost/Infection/DummyInfection.h"
#include "WithinHost/Infection/EmpiricalInfection.h"
#include "WithinHost/Infection/MolineauxInfection.h"
#include "WithinHost/Infection/PennyInfection.h"
#include "WithinHost/Pathogenesis/PathogenesisModel.h"
#include "WithinHost/Diagnostic.h"
#include "inputData.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
//using namespace std;

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

void WHFalciparum::init() {
    sigma_i=sqrt(InputData.getParameter(Params::SIGMA_I_SQ));
    immPenalty_22=1-exp(InputData.getParameter(Params::IMMUNITY_PENALTY));
    immEffectorRemain=exp(-InputData.getParameter(Params::IMMUNE_EFFECTOR_DECAY));
    asexImmRemain=exp(-InputData.getParameter(Params::ASEXUAL_IMMUNITY_DECAY));
    
    _ylagLen = TimeStep::intervalsPer5Days.asInt() * 4;
    
    //NOTE: should also call cleanup() on the PathogenesisModel, but it only frees memory which the OS does anyway
    Pathogenesis::PathogenesisModel::init();
    
    /*
    The detection limit (in parasites/ul) is currently the same for PCR and for microscopy
    TODO: in fact the detection limit in Garki should be the same as the PCR detection limit
    The density bias allows the detection limit for microscopy to be higher for other sites
    */
    double densitybias;
    if (util::ModelOptions::option (util::GARKI_DENSITY_BIAS)) {
        densitybias=InputData.getParameter(Params::DENSITY_BIAS_GARKI);
    } else {
        if ((InputData().getAnalysisNo() >= 22) && (InputData().getAnalysisNo() <= 30)) {
            cerr << "Warning: these analysis numbers used to mean use Garki density bias. If you do want to use this, specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong." << endl;
        }
        densitybias=InputData.getParameter(Params::DENSITY_BIAS_NON_GARKI);
    }
    double detectionLimit=InputData().getMonitoring().getSurveys().getDetectionLimit()*densitybias;
    Diagnostic::default_.setDeterministic( detectionLimit );
    
    Infection::init();
    
    if (util::ModelOptions::option (util::DUMMY_WITHIN_HOST_MODEL)) {
        DummyInfection::init ();
    } else if (util::ModelOptions::option (util::EMPIRICAL_WITHIN_HOST_MODEL)) {
        EmpiricalInfection::init();     // 1-day timestep check
    } else if (util::ModelOptions::option (util::MOLINEAUX_WITHIN_HOST_MODEL)) {
        MolineauxInfection::init();
    } else if (util::ModelOptions::option (util::PENNY_WITHIN_HOST_MODEL)) {
        PennyInfection::init();
    } else {
        DescriptiveInfection::init ();  // 5-day timestep check
    }
}


// -----  Non-static  -----

WHFalciparum::WHFalciparum():
    WHInterface(),
    _cumulativeh(0.0), _cumulativeY(0.0), _cumulativeYlag(0.0),
    totalDensity(0.0), timeStepMaxDensity(0.0)
{
    _innateImmSurvFact = exp(-random::gauss(0, sigma_i));
    
    _ylag.assign (_ylagLen, 0.0);
}
void WHFalciparum::setComorbidityFactor(double factor)
{
    pathogenesisModel = auto_ptr<Pathogenesis::PathogenesisModel>(
        Pathogenesis::PathogenesisModel::createPathogenesisModel(factor) );
}

WHFalciparum::~WHFalciparum()
{
}

double WHFalciparum::probTransmissionToMosquito( TimeStep ageTimeSteps, double tbvEfficacy ) const{
    /* This model (often referred to as the gametocyte model) was designed for
    5-day timesteps. We use the same model (sampling 10, 15 and 20 days ago)
    for 1-day timesteps to avoid having to design and analyse a new model.
    Description: AJTMH pp.32-33 */
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
    double probTransmissionToMosquito = transmit * (1.0 - tbvEfficacy );
    util::streamValidate( probTransmissionToMosquito );
    return probTransmissionToMosquito;
}

bool WHFalciparum::diagnosticDefault() const{
    return Diagnostic::default_.isPositive( totalDensity );
}
bool WHFalciparum::diagnosticMDA() const{
    return Diagnostic::mda.isPositive( totalDensity );
}
Pathogenesis::State WHFalciparum::determineMorbidity(double ageYears){
    Pathogenesis::State state =
            pathogenesisModel->determineState( ageYears, timeStepMaxDensity, totalDensity );
    
    if( (state & Pathogenesis::MALARIA) && util::ModelOptions::option( util::PENALISATION_EPISODES ) ){
        // This does immunity penalisation:
        _cumulativeY = _cumulativeYlag - immPenalty_22*(_cumulativeY-_cumulativeYlag);
        if (_cumulativeY < 0) {
            _cumulativeY=0.0;
        }
    }
    
    return state;
}


// -----  immunity  -----

void WHFalciparum::updateImmuneStatus() {
    if (immEffectorRemain < 1) {
        _cumulativeh*=immEffectorRemain;
        _cumulativeY*=immEffectorRemain;
    }
    if (asexImmRemain < 1) {
        _cumulativeh*=asexImmRemain/
                      (1+(_cumulativeh*(1-asexImmRemain)/Infection::cumulativeHstar));
        _cumulativeY*=asexImmRemain/
                      (1+(_cumulativeY*(1-asexImmRemain)/Infection::cumulativeYstar));
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
