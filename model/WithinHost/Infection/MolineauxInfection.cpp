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
// All the values of q^1... q^v are stored in this array.
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

// Blood volume per kg can be found, for example, here:
// http://hypertextbook.com/facts/1998/LanNaLee.shtml
// Futher searches give values between around 65 and 85 ml/kg.
const double blood_vol_per_kg = 7e4;    // μl/kg; (equals 70 ml/kg)
// PRBC = parasited red blood cell
const double initial_dens = 0.1;        // initial parasite density in PRBC/μl
const double elim_parasites = 50;       // minimum number PRBC to avoid elimination of infection and variant

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
            //TODO: review this variant (or discard)?
            first_local_maximum_gamma = true;
            mean_shape_first_local_max = pow(mean_shape_first_local_max,2) / pow(sd_scale_first_local_max,2);
            sd_scale_first_local_max = pow(sd_scale_first_local_max,2) / mean_shape_first_local_max;
        } else {
            first_local_maximum_gamma = false;
        }
        
        // with gamma distribution shape and scale parameters has to be recalculated
        if(util::ModelOptions::option (util::MEAN_DURATION_GAMMA)) {
            //TODO: review this variant (or discard)?
            mean_duration_gamma = true;
            mean_shape_diff_pos_days = pow(mean_shape_diff_pos_days,2) / pow(sd_scale_diff_pos_days,2);
            sd_scale_diff_pos_days = pow(sd_scale_diff_pos_days,2) / mean_shape_diff_pos_days;
        } else {
            mean_duration_gamma = false;
        }
    }
    
   for( size_t i = 0; i < v; i++ ){
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
        P(0.0), P1(0.0), P2(0.0), variantSpecificSummation(0.0)
{
    for (size_t tau=0; tau<taus; tau++){
        laggedP[tau] =  0.0;
    }
}

// ———  MolineauxInfection: density updates  ———

bool MolineauxInfection::updateDensity( double survivalFactor, SimTime bsAge, double body_mass ){
    // bsAge : age of blood stage; 0 implies initial density (0.1; time t=0 in
    // paper, t=1 in MP's Matlab code), age 2 days is after first update step
    // survivalFactor : probabilty of merozoites surviving drugs, inter-infection immunity and vaccines
    
    double blood_volume = blood_vol_per_kg * body_mass;
    double elim_dens = elim_parasites / blood_volume;   // 50 / 5e6 = 5e-5
    
    if (bsAge == sim::zero()){
        // The first variant starts with a pre-set density (regardless of blood
        // volume; this is an assumption by DH; paper assumes fixed volume)
        variants.resize(1);
        variants[0].P = initial_dens;
        m_density = initial_dens;
    }else{
        double sum = 0.0;
        for( size_t i=0; i<variants.size(); i++ ){
            double newP = survivalFactor * variants[i].P1;
            variants[i].P = static_cast<float>(newP);
            variants[i].P1 = static_cast<float>(survivalFactor * variants[i].P2);
            sum += newP;
        }
        m_density = sum;
    }
    
    // Integration with inter-infection immunity model
    m_cumulativeExposureJ += m_density;
    
    if( m_density <= elim_dens ){
        return true;    // infection goes extinct
    }
    
    // If the infection hasn't gone extinct age is even, then update
    // growthRateMultiplier for the next two steps:
    if( mod_nn(bsAge.inDays(), 2) == 0 ){
        updateGrowthRateMultiplier( bsAge.inDays(), elim_dens );
    }
    return false;
}


void MolineauxInfection::updateGrowthRateMultiplier( int ageDays, double elim_dens ){
    // The immune responses are represented by the variables
    // Sc (probability that a parasite escapes control by innate and variant-transcending immune response)
    // Sm (                        "                      acquired and variant-transcending immune response)
    // S[i] (                      "                      acquired and variant-specific immune response)
    
    // ———  1. innate immunity (equation 5)  ———
    BOOST_STATIC_ASSERT( kappa_c == 3 );        // allows us to optimise pow(..., kappa_c) to multiplication
    const double base = m_density/Pc_star;
    const double Sc = 1.0 / (1.0 + base*base*base);
    
    // ———  2. variant-transcending immune response (equations 7, 8)  ———
    // 2.a) Update the sum in (7), based on previous result
    const size_t tau = mod_nn(ageDays / 2, taus);        // 8 days ago has same index as today
    // NOTE: rho == 0 so we can optimise this code:
    //variantTranscendingSummation = variantTranscendingSummation * exp(-2.0*rho) + laggedPc[index];
    variantTranscendingSummation = variantTranscendingSummation + laggedPc[tau];
    
    // 2.b) Update laggedPc (equation 8, but only stored for the last tau steps)
    //We could use min here, but it seems that min has problems with static const double C
    laggedPc[tau] = static_cast<float>(m_density < C ? m_density : C);
    
    // 2.c) calculate S_m(t) (equation 7)
    BOOST_STATIC_ASSERT( kappa_m == 1 );        // allows us to optimise out pow(..., kappa_m):
    //double Sm = (1.0 - beta) / (1.0 + pow(variantTranscendingSummation / Pm_star, kappa_m)) + beta
    const double Sm = (1.0 - beta) / (1.0 + variantTranscendingSummation / Pm_star) + beta;
    
    // ———  3. variant-specific immune response (equation 6)  ———
    double S[v];        // calculate value for each variant
    double sum_qj_Sj=0.0;       // simultaneously calculation summataion in equation 4
    
    for (size_t i=0; i<v; i++){
        if( i < variants.size() ){
            // 3.a) Update the sum in (6) based on the last step's value
            //note: sigma_decay = exp(-2*sigma)
            variants[i].variantSpecificSummation = static_cast<float>(
                variants[i].variantSpecificSummation * sigma_decay + variants[i].laggedP[tau]);
            // 3.b) update history of density (P_i(t))
            variants[i].laggedP[tau] = variants[i].P;
            
            // 3.c) calculate S_i(t) (equation 6)
            BOOST_STATIC_ASSERT( kappa_v == 3 );        // again, optimise pow to multiplication
            const double base = variants[i].variantSpecificSummation * inv_Pstar_v;
            S[i] = 1.0 / (1.0 + base*base*base);        // eqn 6, given κ_v = 3
        }else{
            S[i] = 1.0; // eqn 6 for the case when P_i(τ) = 0 for τ ≤ t - δ_m
        }
        sum_qj_Sj += qPow[i] * S[i];    // sum in eqn 4
    }
    
    // ———  4. Variant densities, equations 1, 2 and 4  ———
    for (size_t i=0;i<v;i++)
    {
        // 4.a) Calculate p_i, variant selection probability (eqn 4)
        double p_i = 0.0;
        if( S[i] >= 0.1 ){
            //note: qPow[i] = pow(q, i+1)
            p_i = qPow[i] * S[i] / sum_qj_Sj;
        }
        
        // 4.b) calculate P_i'(t+2) [eqn 1] then P_i(t+2) [eqn 2]
        // This is the growth rate after taking immune effect into account:
        double growth_factor = m[i] * S[i] * Sc * Sm;   // part of eqn 1
        if( i < variants.size() ){
            // P_prime: the variant's density at time t+2 (eqn 1)
            double P_prime = ( (1.0 - sProb) * variants[i].P + sProb * p_i * m_density ) * growth_factor;
            
            if( P_prime < elim_dens ) P_prime = 0.0;    // eqn 2
            
            variants[i].P1 = static_cast<float>(sqrt(variants[i].P * P_prime));
            variants[i].P2 = static_cast<float>(P_prime);
        }else{
            // In this case P_i(τ) = 0 for all τ ≤ t, and we haven't allocated storage.
            
            // P_prime: the variant's density at time t+2 (eqn 1 in paper)
            double P_prime = ( sProb * p_i * m_density ) * growth_factor;
            
            // Molineaux paper equation 2
            if( P_prime >= elim_dens ){    // [if not, P_i(t+2) = 0]
                // express a new variant at time t+2:
                variants.resize( i+1 ); // allocate (potentially for multiple variants)
                variants[i].P2 = static_cast<float>(P_prime);
            }
        }
    }
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
        P & stream;
        P1 & stream;
        P2 & stream;
        variantSpecificSummation & stream;
        for(size_t i = 0; i < taus; ++i){
            laggedP[i] & stream;
        }
    }
    // else: all members are zero-initialised by the constructor, so don't do anything
}

void MolineauxInfection::Variant::operator& (ostream& stream) {
    bool nonZero =
            P != 0.0 ||
            P1 != 0.0 ||
            P2 != 0.0 ||
            variantSpecificSummation != 0.0;

    nonZero & stream;
    if( nonZero ){
        P & stream;
        P1 & stream;
        P2 & stream;
        variantSpecificSummation & stream;
        for(size_t i = 0; i < taus; ++i){
            laggedP[i] & stream;
        }
    }
}

}
}
