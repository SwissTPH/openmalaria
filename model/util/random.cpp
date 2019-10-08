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

/* This module contains the random-number generator and distributions wrapper.
 *
 * Currently both the GSL and boost generators are implemented. The
 * distributions all come from the GSL library so far.
 * 
 * Using the boost generator appears (in rough tests) to be slightly
 * slower, which is understandable since the GSL distributions must then use a
 * wrapper around the boost generator.
 * 
 * Note: using boost distributions elsewhere could ideally be implemented a
 * little differently, since the distribution objects could in many cases last
 * the length of the program rather than be created on each use.
 */

// Define to use boost distributions. Unfortunately these are not necessarily
// value-stable or portable.
// #define OM_RANDOM_USE_BOOST_DIST

#include "util/random.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "Global.h"

#ifdef OM_RANDOM_USE_BOOST_DIST
#include <boost/random/normal_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/weibull_distribution.hpp>
#endif

#if !defined OM_RANDOM_USE_BOOST_DIST
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#endif

#include <cmath>
#include <sstream>

// Note: since we're using both gsl and boost files, we should be careful to
// avoid name conflicts. So probably don't use "using namespace boost;".


namespace OM { namespace util {

RNG global_RNG;

// Support functions
long unsigned int sample_ulong (void *ptr) {
    pcg32 *rng = reinterpret_cast<pcg32*>(ptr);
    BOOST_STATIC_ASSERT (sizeof(uint32_t) <= sizeof(long unsigned int));
    return static_cast<long unsigned int> ((*rng)());
}
double sample_double01 (void *ptr) {
    pcg32 *rng = reinterpret_cast<pcg32*>(ptr);
    boost::random::uniform_01<pcg32&> dist (*rng);
    return dist ();
}

gsl_rng_type make_gsl_rng_type(pcg32& rng) {
    return {
        "PCG32",			// name
        rng.max(),		// max value
        rng.min(),		// min value
        0,				// size of state; not used here
        nullptr,			// re-seed function; don't use
        &sample_ulong,
        &sample_double01
    };
}


RNG::RNG() {
    m_gsl_type = make_gsl_rng_type(m_rng);
    m_gsl_gen.type = &m_gsl_type;
    m_gsl_gen.state = reinterpret_cast<void*>(&m_rng);
}

RNG::RNG(RNG&& other) {
    m_rng = other.m_rng;
    m_gsl_type = other.m_gsl_type;
    // The following point into instance of this RNG:
    m_gsl_gen.type = &m_gsl_type;
    m_gsl_gen.state = reinterpret_cast<void*>(&m_rng);
}
void RNG::operator=(RNG&& other) {
    m_rng = other.m_rng;
    m_gsl_type = other.m_gsl_type;
    // The following point into instance of this RNG:
    m_gsl_gen.type = &m_gsl_type;
    m_gsl_gen.state = reinterpret_cast<void*>(&m_rng);
}

void RNG::seed (uint64_t seed) {
    m_rng.seed(seed);
}

// Incomplete checkpointing: assumes seed(..) has been called as before
void RNG::checkpoint (ostream& stream) {
    m_rng.binary_checkpoint(stream);
}
void RNG::checkpoint(istream& stream) {
    m_rng.binary_checkpoint(stream);
}



// -----  random number generation  -----

double RNG::gauss (double mean, double std){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::normal_distribution<> dist (mean, std);
    double result = dist(m_rng);
# else
    double result = gsl_ran_gaussian(&m_gsl_gen,std)+mean;
# endif
//     util::streamValidate(result);
    return result;
}

double RNG::gamma (double a, double b){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::gamma_distribution<> dist (a, b);
    double result = dist(m_rng);
# else
    double result = gsl_ran_gamma(&m_gsl_gen, a, b);
# endif
//     util::streamValidate(result);
    return result;
}

double RNG::log_normal (double meanlog, double stdlog){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::lognormal_distribution<> dist (meanlog, stdlog);
    double result = dist (m_rng);
# else
    double result = gsl_ran_lognormal (&m_gsl_gen, meanlog, stdlog);
# endif
//     util::streamValidate(result);
    return result;
}

double RNG::max_multi_log_normal (double start, int n, double meanlog, double stdlog){
    double result = start;
# ifdef OM_RANDOM_USE_BOOST_DIST
    // We don't have a CDF available, so take log-normal samples
    for (int i = 0; i < n; ++i) {
        double sample = log_normal(meanlog, stdlog);
        if (sample > result) result = sample;
    }
# else
    // Used for performance reasons. Calling GSL's log_normal 5 times is 50% slower.
    // For random variables X1, .., Xn with identical distribution X and
    // Mn = max(X1, .., Xn), the CDF:
    //      F_Mn(x) = P(Mn ≤ x)
    //              = P(X1 ≤ x) P(X2 ≤ x) .. P(Xn ≤ x)
    //              = F_X(x) ^ n
    // Thus for u = F_Mn(x) = F_X(x) ^ n, u^(1/n) = F_X(x).
    double normp = pow( uniform_01(), 1.0 / n );
    double zval = gsl_cdf_ugaussian_Pinv (normp);
    double multi_sample = exp(meanlog + stdlog * zval);
    result = std::max(result, multi_sample);
# endif
//     util::streamValidate(result);
    return result;
}

double RNG::beta (double a, double b){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::beta_distribution<> dist (a, b);
    double result = dist(m_rng);
# else
    double result = gsl_ran_beta (&m_gsl_gen,a,b);
# endif
//     util::streamValidate(result);
    return result;
}

int RNG::poisson(double lambda){
    if( !(boost::math::isfinite)(lambda) ){
	//This would lead to an inifinite loop
	throw TRACED_EXCEPTION( "lambda is inf", Error::InfLambda );
    }
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::poisson_distribution<> dist (lambda);
    int result = dist(m_rng);
# else
    int result = gsl_ran_poisson (&m_gsl_gen, lambda);
# endif
//     util::streamValidate(result);
    return result;
}

bool RNG::bernoulli(double prob){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::bernoulli_distribution<> dist (prob);
    bool result = dist(m_rng);
# else
    assert( (boost::math::isfinite)(prob) );
    // return true iff our variate is less than the probability
    bool result = uniform_01() < prob;
# endif
//     util::streamValidate(result);
    return result;
}

double RNG::weibull(double lambda, double k){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::weibull_distribution<> dist (k, lambda);
    return dist(m_rng);
# else
    return gsl_ran_weibull( &m_gsl_gen, lambda, k );
# endif
}

} }
