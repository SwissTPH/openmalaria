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
#include "WithinHost/Genotypes.h"
#include "mon/reporting.h"
#include "util/random.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/checkpoint_containers.h"
#include "util/timeConversions.h"
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
int WHFalciparum::y_lag_len = 0;

// -----  static functions  -----

void WHFalciparum::init( const OM::Parameters& parameters, const scnXml::Model& model ) {
    sigma_i=sqrt(parameters[Parameters::SIGMA_I_SQ]);
    immPenalty_22=1-exp(parameters[Parameters::IMMUNITY_PENALTY]);
    immEffectorRemain=exp(-parameters[Parameters::IMMUNE_EFFECTOR_DECAY]);
    asexImmRemain=exp(-parameters[Parameters::ASEXUAL_IMMUNITY_DECAY]);
    
    y_lag_len = SimTime::daysToSteps(20);
    
    //NOTE: should also call cleanup() on the PathogenesisModel, but it only frees memory which the OS does anyway
    Pathogenesis::PathogenesisModel::init( parameters, model.getClinical(), false );
    
    try{
        //NOTE: if XSD is changed, this should not have a default unit
        SimTime latentP = UnitParse::readShortDuration(
            model.getParameters().getLatentp(), UnitParse::STEPS );
        Infection::init( parameters, latentP );
    }catch( const util::format_error& e ){
        throw util::xml_scenario_error( string("model/parameters/latentP: ").append(e.message()) );
    }
}


// -----  Non-static  -----

WHFalciparum::WHFalciparum( double comorbidityFactor ):
    WHInterface(),
    m_cumulative_h(0.0), m_cumulative_Y(0.0), m_cumulative_Y_lag(0.0),
    totalDensity(0.0), timeStepMaxDensity(0.0),
    pathogenesisModel( Pathogenesis::PathogenesisModel::createPathogenesisModel( comorbidityFactor ) )
{
    // NOTE: negating a Gaussian sample with mean 0 is pointless — except that
    // the individual samples change. In any case the overhead is negligible.
    //FIXME: Should this be allowed to be greater than 1?
    // Oldest code on GoogleCode: _innateImmunity=(double)(W_GAUSS((0), (sigma_i)));
    _innateImmSurvFact = exp(-random::gauss(sigma_i));
    
    m_y_lag.assign(y_lag_len, Genotypes::N(), 0.0);
}

WHFalciparum::~WHFalciparum()
{
}

// Infectiousness parameters: see AJTMH p.33; tau=1/sigmag**2 
const double PTM_beta1=1.0;
const double PTM_beta2=0.46;
const double PTM_beta3=0.17;
const double PTM_tau= 0.066;
const double PTM_tau_prime = 1.0 / sqrt(1.0 / PTM_tau);
const double PTM_mu= -8.1;

double WHFalciparum::probTransmissionToMosquito( double tbvFactor, double *sumX ) const{
    // This model (often referred to as the gametocyte model) was designed for
    // 5-day time steps. We use the same model (sampling 10, 15 and 20 days
    // ago) for 1-day time steps to avoid having to design and analyse a new
    // model. Description: AJTMH pp.32-33 and p9.
    
    // Note: we don't allow for gametocydal treatments (e.g. Primaquine).
    
    // Take weighted sum of total asexual blood stage density 10, 15 and 20 days
    // before. Add y_lag_len to index to ensure positive.
    const int i10 = (sim::ts0() - SimTime::fromDays(10) + SimTime::oneTS()).inSteps() + y_lag_len;
    const int i5d = SimTime::daysToSteps(5);
    const int i10d = 2 * i5d;
    // Sum lagged densities across genotypes:
    double y10 = 0.0, y15 = 0.0, y20 = 0.0;
    for( size_t genotype = 0; genotype < Genotypes::N(); ++genotype ){
        y10 += m_y_lag.at(mod_nn(i10, y_lag_len), genotype);
        y15 += m_y_lag.at(mod_nn(i10 - i5d, y_lag_len), genotype);
        y20 += m_y_lag.at(mod_nn(i10 - i10d, y_lag_len), genotype);
    }
    // Weighted sum:
    const double x = PTM_beta1 * y10 + PTM_beta2 * y15 + PTM_beta3 * y20;
    if( sumX != 0 ) *sumX = 1.0 / x;    // copy to sumX, if set
    if( x < 0.001 ) return 0.0; // cut off for uninfectious humans
    
    // Get a zval, convert to equivalent Normal sample:
    const double zval = (log(x) + PTM_mu) * PTM_tau_prime;
    const double pone = gsl_cdf_ugaussian_P(zval);
    double pTransmit = pone*pone;
    // pTransmit has to be between 0 and 1:
    pTransmit=std::max(pTransmit, 0.0);
    pTransmit=std::min(pTransmit, 1.0);
    
    // Include here the effect of transmission-blocking vaccination:
    pTransmit *= tbvFactor;
    util::streamValidate( pTransmit );
    return pTransmit;
}
double WHFalciparum::pTransGenotype(double pTrans, double sumX, size_t genotype)
{
    assert( pTrans > 0.0 );
    assert( (boost::math::isfinite)(sumX) );
    
    // This is an extension of the original model.
    //NOTE: it is an approximation since it ignores the possibility of
    // simultaneously infecting a mosquito with multiple genotypes.
    
    // Take weighted sum of total asexual blood stage density 10, 15 and 20 days
    // before. Add y_lag_len to index to ensure positive.
    const int i10 = (sim::ts0() - SimTime::fromDays(10) + SimTime::oneTS()).inSteps() + y_lag_len;
    const int i5d = SimTime::daysToSteps(5);
    const int i10d = 2 * i5d;
    const double x =
        PTM_beta1 * m_y_lag.at(mod_nn(i10, y_lag_len), genotype) +
        PTM_beta2 * m_y_lag.at(mod_nn(i10 - i5d, y_lag_len), genotype) +
        PTM_beta3 * m_y_lag.at(mod_nn(i10 - i10d, y_lag_len), genotype);
    
    return pTrans * x * sumX;
}

bool WHFalciparum::diagnosticResult( const Diagnostic& diagnostic ) const{
    return diagnostic.isPositive( totalDensity );
}

void WHFalciparum::treatment( Host::Human& human, TreatmentId treatId ){
    const Treatments& treat = Treatments::select( treatId );
    treatSimple( human, treat.liverEffect(), treat.bloodEffect() );
    
    // triggered intervention deployments:
    treat.deploy( human,
                  mon::Deploy::TREAT,
                  interventions::VaccineLimits(/*default initialise: no limits*/) );
}
bool WHFalciparum::treatSimple( const Host::Human& human, SimTime timeLiver, SimTime timeBlood ){
    if( timeLiver != SimTime::zero() ){
        if( timeLiver < SimTime::zero() )
            clearInfections( Treatments::LIVER );
        else
            treatExpiryLiver = max( treatExpiryLiver, sim::nowOrTs1() + timeLiver );
        mon::reportEventMHI( mon::MHT_LS_TREATMENTS, human, 1 );
    }
    if( timeBlood != SimTime::zero() ){
        if( timeBlood < SimTime::zero() )
            clearInfections( Treatments::BLOOD );
        else
            treatExpiryBlood = max( treatExpiryBlood, sim::nowOrTs1() + timeBlood );
        return true;    // blood stage treatment
    }
    return false;    // no blood stage treatment
}

Pathogenesis::StatePair WHFalciparum::determineMorbidity( Host::Human& human,
        double ageYears, bool isDoomed )
{
    Pathogenesis::StatePair result = pathogenesisModel->determineState( human,
            ageYears, timeStepMaxDensity, totalDensity, isDoomed );
    
    /* Note: this model can easily be re-enabled, but is not used and not considered to be a good model.
    if( (result.state & Pathogenesis::MALARIA) && util::ModelOptions::option( util::PENALISATION_EPISODES ) ){
        // This does immunity penalisation:
        m_cumulative_Y = m_cumulative_Y_lag - immPenalty_22*(m_cumulative_Y-m_cumulative_Y_lag);
        if (m_cumulative_Y < 0) {
            m_cumulative_Y=0.0;
        }
    }*/
    
    return result;
}


// -----  immunity  -----

void WHFalciparum::updateImmuneStatus() {
    if (immEffectorRemain < 1) {
        m_cumulative_h*=immEffectorRemain;
        m_cumulative_Y*=immEffectorRemain;
    }
    if (asexImmRemain < 1) {
        m_cumulative_h*=asexImmRemain/
                      (1+(m_cumulative_h*(1-asexImmRemain) * Infection::invCumulativeHstar));
        m_cumulative_Y*=asexImmRemain/
                      (1+(m_cumulative_Y*(1-asexImmRemain) * Infection::invCumulativeYstar));
    }
    m_cumulative_Y_lag = m_cumulative_Y;
}


// -----  Checkpointing  -----

void WHFalciparum::checkpoint (istream& stream) {
    WHInterface::checkpoint( stream );
    _innateImmSurvFact & stream;
    m_cumulative_h & stream;
    m_cumulative_Y & stream;
    m_cumulative_Y_lag & stream;
    totalDensity & stream;
    timeStepMaxDensity & stream;
    m_y_lag & stream;
    (*pathogenesisModel) & stream;
    treatExpiryLiver & stream;
    treatExpiryBlood & stream;
}
void WHFalciparum::checkpoint (ostream& stream) {
    WHInterface::checkpoint( stream );
    _innateImmSurvFact & stream;
    m_cumulative_h & stream;
    m_cumulative_Y & stream;
    m_cumulative_Y_lag & stream;
    totalDensity & stream;
    timeStepMaxDensity & stream;
    m_y_lag & stream;
    (*pathogenesisModel) & stream;
    treatExpiryLiver & stream;
    treatExpiryBlood & stream;
}

}
}
