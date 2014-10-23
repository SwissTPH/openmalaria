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

#include "WithinHost/Infection/MolineauxInfection.h"
#include "WithinHost/CommonWithinHost.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/checkpoint_containers.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <boost/static_assert.hpp>


namespace OM {
namespace WithinHost {

using namespace OM::util;

// ———  model constants set from input parameters  ———

// depends on lognormal/gamma distribution for first_local_max
static bool first_local_maximum_gamma = false;
// mean/shape and sd/scale of the first local maximum density for lognormal/gamma distribution
static double mean_shape_first_local_max = std::numeric_limits<double>::signaling_NaN();
static double sd_scale_first_local_max = std::numeric_limits<double>::signaling_NaN();

// depends on between lognormal/gamma distribution for mean_duration diff_pos_days
static bool mean_duration_gamma = false;
// mean/shape and sd/scale of the difference between last positive and first positive days for lognormal/gamma distribution
static double mean_shape_diff_pos_days = std::numeric_limits<double>::signaling_NaN();
static double sd_scale_diff_pos_days = std::numeric_limits<double>::signaling_NaN();

// boolean choosing between gamma and lognormal distribution for equation 11
static bool multi_factor_gamma = false;

// pairwise sample of case-specific P* parameters
static bool pairwise_PStar_sample = false;

// q^(i+1) array
// All the values of q^1... q^50 are stored in this array.
// This avoids the recalculation of those values every second time step. */
static double qPow[MolineauxInfection::v];

// ———  hard-coded model constants  ———

// sigma, rho: decay parameters, per day, of the acquired variant-specific and variant-transcending immune responses
const double sigma = 0.02;
const double sigma_decay = exp(-2.0*sigma);     //FIXME: should not be factor of two here, unless sigma is decay per day!
// const double rho = 0.0;      // NOTE: rho is not used; if setting this, also change the code below (search rho)
// beta: Minimum value of the probability that a parasite escapes control by the acquired and variant-transcending immune response
const double beta = 0.01;
// sProb: the fraction of parasites switching among variants per two-day cycle
const double sProb = 0.02;
// q: Parameter of the geometric distribution of switching probabilities
const double q = 0.3;

// mu_m, sigma_m: Mean and standard deviation to use for the normal distribution setting the variant specific multiplication factor.
const double mu_m = 16.0;
const double sigma_m = 10.4;
// constants for sampling multiplication factors from a Gamma distribution instead of normal (above)
const double shape_m = 2.4;
const double scale_m = 6.8;

// k_c,k_m: constants allowing calculation of Pc_star and Pm_star from host-specific data
const double k_c = 0.2;
const double k_m = 0.04;
// Pstar_v: critical density of a variant, common to all variants
const double inv_Pstar_v = 1.0 / 30.0;
// kappa_c, kappa_m, kappa_v: Stiffness parameters for saturation of immune responses
const int kappa_c = 3;
const int kappa_m = 1;
const int kappa_v = 3;
// C: Maximum daily antigenic stimulus, per mul, of the acquired variant-transcending immune response
const double C = 1.0;

// Density of first variant at start of infection
const float initial_P = 0.1f;

/* Case-specific parameters.
 * 
 * For each of the 35 Malaria Therapy patients, this is:
 *  0) the duration of the infection as last day with positive density minus first positive
 *  1) the parasite density in parasites/microlitre at the first local maximum
 * 
 * Indexing is 2*patient + param where 0 <= patient < 35 and 0 <= param < 2.
 */
const double caseSpecificData[] = {
        216, 18600, // G131
        198, 13080, // G140
        206, 45720, // G141
        366, 23760, // G142
        230, 60840, // G143
        172, 6000, // G146
        100, 2340, // G147
        236, 31440, // G159
        236, 453600, // G161
        120, 4240, // G173
        176, 195840, // G174
        178, 60120, // G178
        36, 8720, // G184
        44, 8000, // G193
        242, 395280, // G200
        70, 28320, // G210
        292, 200160, // G212
        248, 59320, // G217
        98, 66480, // G23
        176, 61200, // G230
        234, 169920, // G240
        226, 46800, // G264
        270, 19260, // G265
        278, 86040, // G279
        212, 110160, // G290
        264, 43200, // G385
        364, 133920, // G407
        184, 222480, // G408
        160, 21420, // G414
        220, 74160, // G416
        132, 210960, // G423
        176, 89280, // G439
        208, 105840, // G445
        330, 21600, // G457
        404, 156240 // G48
};

// ———  static (non-class-member) code  ———

CommonInfection* createMolineauxInfection (uint32_t protID) {
    return new MolineauxInfection (protID);
}

CommonInfection* checkpointedMolineauxInfection (istream& stream) {
    return new MolineauxInfection (stream);
}

void MolineauxInfection::init( const Parameters& parameters ){
    CommonWithinHost::createInfection = &createMolineauxInfection;
    CommonWithinHost::checkpointedInfection = &checkpointedMolineauxInfection;
    
    multi_factor_gamma = util::ModelOptions::option( util::PARASITE_REPLICATION_GAMMA );
    
    pairwise_PStar_sample = util::ModelOptions::option( util::MOLINEAUX_PAIRWISE_SAMPLE );
    if( !pairwise_PStar_sample ){
        mean_shape_first_local_max = parameters[Parameters::MEAN_LOCAL_MAX_DENSITY];
        sd_scale_first_local_max = parameters[Parameters::SD_LOCAL_MAX_DENSITY];
        
        mean_shape_diff_pos_days = parameters[Parameters::MEAN_DIFF_POS_DAYS];
        sd_scale_diff_pos_days = parameters[Parameters::SD_DIFF_POS_DAYS];
        
        
        // with gamma distribution shape and scale parameters has to be recalculated 
        if (util::ModelOptions::option (util::FIRST_LOCAL_MAXIMUM_GAMMA)) {
            first_local_maximum_gamma = true;
            mean_shape_first_local_max = pow(mean_shape_first_local_max,2) / pow(sd_scale_first_local_max,2);
            sd_scale_first_local_max = pow(sd_scale_first_local_max,2) / mean_shape_first_local_max;
        } else {
            first_local_maximum_gamma = false;
        }
        
        // with gamma distribution shape and scale parameters has to be recalculated
        if(util::ModelOptions::option (util::MEAN_DURATION_GAMMA)) {
            mean_duration_gamma = true;
            mean_shape_diff_pos_days = pow(mean_shape_diff_pos_days,2) / pow(sd_scale_diff_pos_days,2);
            sd_scale_diff_pos_days = pow(sd_scale_diff_pos_days,2) / mean_shape_diff_pos_days;
        } else {
            mean_duration_gamma = false;
        }
    }
    
   for(int i=0;i<50;i++){
       qPow[i] = pow(q, static_cast<double>(i+1));
   }
}

// ———  MolineauxInfection: initialisation  ———

MolineauxInfection::MolineauxInfection(uint32_t protID):
        CommonInfection(protID)
{
    for( size_t i=0;i<v; i++ ){
        // Molineaux paper, equation 11
        if( multi_factor_gamma ){
            do{
                m[i] = static_cast<float>(random::gamma(shape_m,scale_m));
            }while( m[i]<1.0 );
        }else{
            do{
                m[i] = static_cast<float>(random::gauss(mu_m, sigma_m));
            }while( m[i]<1.0 );
        }
    }
    
    for( size_t tau=0; tau<taus; tau++ ){
        laggedPc[tau] = 0.0;
    }
    
    variantTranscendingSummation = 0.0;
    
    if( pairwise_PStar_sample ){
        int patient = util::random::uniform( 35 );
        Pc_star = k_c * caseSpecificData[2 * patient + 1];
        Pm_star = k_m * caseSpecificData[2 * patient + 0];
    }else{
        // Sampling of the first local maxima:
        if( first_local_maximum_gamma ){
            Pc_star = static_cast<float>( k_c * pow(10.0, random::gamma(
                mean_shape_first_local_max, sd_scale_first_local_max)) );
        } else {
            Pc_star = static_cast<float>( k_c * pow(10.0, random::gauss(
                mean_shape_first_local_max, sd_scale_first_local_max)) );
        }
        
        // Sampling of duration:
        if( mean_duration_gamma ) {
            Pm_star = static_cast<float>( k_m * pow(10.0, random::gamma(
                mean_shape_diff_pos_days,sd_scale_diff_pos_days)) );
        } else {
            Pm_star = static_cast<float>( k_m * pow(10.0, random::gauss(
                mean_shape_diff_pos_days,sd_scale_diff_pos_days)) );
        }
    }
}

MolineauxInfection::Variant::Variant () :
        growthRate(0.0), P(0.0), variantSpecificSummation(0.0), initP(0.0)
{
    for (size_t tau=0; tau<taus; tau++){
        laggedP[tau] =  0.0;
    }
}

// ———  MolineauxInfection: density updates  ———

bool MolineauxInfection::updateDensity( double survivalFactor, SimTime bsAge ){
    // bsAge : age of blood stage; 0 implies initial density (0.1; time t=0 in
    // paper, t=1 in MP's Matlab code), age 2 days is after first update step
    // survivalFactor : probabilty of merozoites surviving drugs, inter-infection immunity and vaccines
    
    if (bsAge == sim::zero()){
        // The first variant starts with density 0.1 while other variants start at zero
        variants.resize(1);
        variants[0].setP( initial_P );
        m_density = initial_P;
    }else{
        double sum = 0.0;
        for( size_t i=0; i<variants.size(); i++ ){
            sum += variants[i].updateDensity( survivalFactor, bsAge.inDays() );
        }
        m_density = sum;
    }
    
    // Integration with inter-infection immunity model
    m_cumulativeExposureJ += m_density;
    
    if( m_density <= 1.0e-5 ){
        return true;    // infection goes extinct
    }
    
    // If the infection hasn't gone extinct age is even, then update
    // growthRateMultiplier for the next two steps:
    if( mod_nn(bsAge.inDays(), 2) == 0 ){
        updateGrowthRateMultiplier(bsAge.inDays());
    }
    return false;
}

double MolineauxInfection::Variant::updateDensity (double survivalFactor, int ageDays) {
    // growthRate:
    // p(t+1) = p(t) * sqrt(p(t+2)/p(t))
    // p(t+2) = p(t+1) * sqrt(p(t+2)/p(t))
    // p(t+2) = p(t) * sqrt(p(t+2)/p(t))^2...
    P *= growthRate;
    
    // survivalFactor: effects of drugs, immunity and vaccines
    // Note: this gives us an additional immune response (besides Si, Sc, Sm);
    // according to MP this is the best thing to do for now.
    P = static_cast<float>(P * survivalFactor);
    initP = static_cast<float>(initP * survivalFactor);
    
    // if t+2: The new variant is now expressed. For already extinct
    // variants this doesn't matter, since initP = 0 for those variants.
    if( P==0 && mod_nn(ageDays, 2) == 0 ){
        P = initP;
    }
    
    return P;
}


void MolineauxInfection::updateGrowthRateMultiplier(int ageDays) {
    // The immune responses are represented by the variables
    // Sc (probability that a parasite escapes control by innate and variant-transcending immune response)
    // Sm (                        "                      acquired and variant-transcending immune response)
    // S[i] (                      "                      acquired and variant-specific immune response)
    
    BOOST_STATIC_ASSERT( kappa_c == 3 );        // allows us to optimise pow(..., kappa_c) to multiplication
    double base = m_density/Pc_star;
    double Sc = 1.0 / (1.0 + base*base*base);
    
    BOOST_STATIC_ASSERT( kappa_m == 1 );        // allows us to optimise out pow(..., kappa_m):
    //double Sm = (1.0 - beta) / (1.0 + pow(getVariantTranscendingSummation(ageDays) * inv_Pm_star, kappa_m)) + beta
    double Sm = (1.0 - beta) / (1.0 + getVariantTranscendingSummation(ageDays) / Pm_star) + beta;
    double S[v];

    double sum_qj_Sj=0.0;
    
    for (size_t i=0; i<v; i++)
    {
        S[i] = 1.0;
        if( i < variants.size() ){
            BOOST_STATIC_ASSERT( kappa_v == 3 );        // again, optimise pow to multiplication
            double base = variants[i].getVariantSpecificSummation(ageDays) * inv_Pstar_v;
            S[i] = 1.0 / (1.0 + base*base*base);
        }
        sum_qj_Sj+= qPow[i]*S[i];
    }

    for (size_t i=0;i<v;i++)
    {
        // Molineaux paper equation 4
        // p_i: variant selection probability
        double p_i = 0.0;
        if( S[i] >= 0.1 ){
            //note: qPow[i] = pow(q, i+1)
            p_i = qPow[i] * S[i] / sum_qj_Sj;
        }
        
        // This is the growth rate after taking immune effect into account:
        double growth_factor = m[i] * S[i] * Sc * Sm;
        if( i < variants.size() ){
            variants[i].updateGrowthRateMultiplier( p_i*m_density, growth_factor );
        }else{
            //NOTE: this does the same as updateGrowthRateMultiplier, where P == 0
            // The code duplication is due to an optimisation: not storing variants which haven't been expressed.
            
            // P_prime: Variant density at t = t + 2 (eqn 1 in paper)
            double P_prime = ( sProb * p_i * m_density ) * growth_factor;
            
            // Molineaux paper equation 2
            if( P_prime >= 1.0e-5 ){    // [if not, Pi(t+2)) = 0 and we don't add a new variant]
                // express a new variant:
                variants.resize( i+1 );
                variants[i].setNextP( static_cast<float>(P_prime) );
            }
        }
    }
}

void MolineauxInfection::Variant::updateGrowthRateMultiplier( double pd, double growth_factor ){
    // P_prime: Variant density at t = t + 2 (eqn 1 in paper)
    double P_prime = ( (1.0 - sProb) * P + sProb * pd ) * growth_factor;
    
    // Molineaux paper equation 2
    if (P_prime<1.0e-5){
        P_prime = 0.0;
    }
    
    // if P == 0 then that means this variant wasn't expressed yet
    // or is extinct. If this variant is emerging in (t+2) the new variant
    // density is stored in the initP array,  so that we are able to add
    // the survival factor's effect to the emerging variant's density.
    if( P == 0 ){
        initP = static_cast<float>(P_prime);
        growthRate = 0.0;
    }else{
        initP = 0.0;
        growthRate = static_cast<float>(sqrt(P_prime / P));
    }
}

double MolineauxInfection::getVariantTranscendingSummation(int ageDays) {
    //Molineaux paper equation 5
    size_t tau = mod_nn(ageDays / 2, taus);        // 8 days ago has same index as today
    // NOTE: rho == 0 so we can optimise this code:
    //variantTranscendingSummation = variantTranscendingSummation * exp(-2.0*rho) + laggedPc[index];
    variantTranscendingSummation = variantTranscendingSummation + laggedPc[tau];
    
    //Molineaux paper equation 8
    //We could use min here, but it seems that min has problems with static const double C
    laggedPc[tau] = static_cast<float>(m_density < C ? m_density : C);
    
    return variantTranscendingSummation;
}

double MolineauxInfection::Variant::getVariantSpecificSummation(int ageDays) {
    //The effective exposure is computed by adding in the 8-day lagged parasite density (i.e. 4 time steps)
    //and decaying the previous value for the effective exposure with decay parameter 2*sigma (the 2 arises because
    //the time steps are two days and the unit of sigma is per day. (reasoning: rearrangment of Molineaux paper equation 6)
    
    //Molineaux paper equation 6
    size_t tau = mod_nn(ageDays / 2, taus);      // 8 days ago has same index as today
    //note: sigma_decay = exp(-2*sigma)
    variantSpecificSummation = static_cast<float>(variantSpecificSummation * sigma_decay + laggedP[tau]);
    laggedP[tau] = P;
    
    return variantSpecificSummation;
}

// ———  MolineauxInfection: checkpointing  ———

MolineauxInfection::MolineauxInfection (istream& stream) :
        CommonInfection(stream)
{
    variantTranscendingSummation & stream;
    for (size_t i=0;i<v;i++) {
        m[i] & stream;
    }
    variants & stream;
    for (size_t j=0;j<taus;j++){
        laggedPc[j] & stream;
    }
    Pc_star & stream;
    Pm_star & stream;
}

void MolineauxInfection::checkpoint (ostream& stream) {
    CommonInfection::checkpoint (stream);

    variantTranscendingSummation & stream;
    for (size_t i=0;i<v;i++) {
        m[i] & stream;
    }
    variants & stream;
    for (size_t j=0;j<taus;j++){
        laggedPc[j] & stream;
    }
    Pc_star & stream;
    Pm_star & stream;
}

void MolineauxInfection::Variant::operator& (istream& stream) {
    bool nonZero;
    nonZero & stream;
    if( nonZero ){
        growthRate & stream;
        P & stream;
        variantSpecificSummation & stream;
        initP & stream;
        for(size_t i = 0; i < taus; ++i){
            laggedP[i] & stream;
        }
    }
    // else: all members are zero-initialised by the constructor, so don't do anything
}

void MolineauxInfection::Variant::operator& (ostream& stream) {
    bool nonZero =
            growthRate != 0.0 ||
            P != 0.0 ||
            variantSpecificSummation != 0.0 ||
            initP != 0.0;

    nonZero & stream;
    if( nonZero ){
        growthRate & stream;
        P & stream;
        variantSpecificSummation & stream;
        initP & stream;
        for(size_t i = 0; i < taus; ++i){
            laggedP[i] & stream;
        }
    }
}

}
}
