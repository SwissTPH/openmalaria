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

#include <boost/random/uniform_01.hpp>
#include <pcg_random.hpp>

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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

#include <cmath>
#include <sstream>

// Note: since we're using both gsl and boost files, we should be careful to
// avoid name conflicts. So probably don't use "using namespace boost;".


namespace OM { namespace util {

    // I would prefer to use pcg64, but MSVC mysteriously fails
    static pcg32 generator;
    static boost::random::uniform_01<pcg32&> rng_uniform01 (generator);
    
#   if !defined OM_RANDOM_USE_BOOST_DIST
    long unsigned int boost_rng_get (void*) {
	BOOST_STATIC_ASSERT (sizeof(uint32_t) <= sizeof(long unsigned int));
	long unsigned int val = static_cast<long unsigned int> (generator ());
	streamValidate( val );
	return val;
    }
    double boost_rng_get_double_01 (void*) {
	return rng_uniform01 ();
    }
    
    static const gsl_rng_type boost_mt_type = {
	"boost_mt19937",		// name
	generator.max(),	// max value
	generator.min(),	// min value
	0,					// size of state; not used here
	nullptr,				// re-seed function; don't use
	&boost_rng_get,
	&boost_rng_get_double_01
    };
#   endif

#if !defined OM_RANDOM_USE_BOOST_DIST
// This should be created and deleted automatically, taking care of
// allocating and freeing the generator.
struct generator_factory {
    gsl_rng * gsl_generator;
    
    generator_factory () {
	// In this case, I construct a wrapper around boost's generator. The reason for this is
	// that it allows use of distributions from both boost and GSL.
	gsl_generator = new gsl_rng;
	gsl_generator->type = &boost_mt_type;
	gsl_generator->state = nullptr;	// state is stored as static variables
    }
    ~generator_factory () {
	delete gsl_generator;
    }
} rng;
# endif

// -----  set-up, tear-down and checkpointing  -----

void random::seed (uint32_t seed) {
//     util::streamValidate(seed);
    generator.seed (seed);
}

void random::checkpoint (istream& stream, int seedFileNumber) {
    // Don't use OM::util::checkpoint function for loading a stream; checkpoint::validateListSize uses too small a number.
    string str;
    size_t len;
    len & stream;
    str.resize (len);
    stream.read (&str[0], str.length());
    if (!stream || stream.gcount() != streamsize(len))
	throw checkpoint_error ("stream read error string");
    istringstream ss (str);
    ss >> generator;
}

void random::checkpoint (ostream& stream, int seedFileNumber) {
    ostringstream ss;
    ss << generator;
    ss.str() & stream;
}


// -----  random number generation  -----

double random::uniform_01 () {
    double result =
    // GSL and boost versions both do the same (when using boost as the underlying generator):
# ifdef OM_RANDOM_USE_BOOST_DIST
        rng_uniform01 ();
# else
        gsl_rng_uniform (rng.gsl_generator);
# endif
//     util::streamValidate(result);
    return result;
}

double random::gauss (double mean, double std){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::normal_distribution<> dist (mean, std);
    double result = dist(generator);
# else
    double result = gsl_ran_gaussian(rng.gsl_generator,std)+mean;
# endif
//     util::streamValidate(result);
    return result;
}

double random::gamma (double a, double b){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::gamma_distribution<> dist (a, b);
    double result = dist(generator);
# else
    double result = gsl_ran_gamma(rng.gsl_generator, a, b);
# endif
//     util::streamValidate(result);
    return result;
}

double random::log_normal (double meanlog, double stdlog){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::lognormal_distribution<> dist (meanlog, stdlog);
    double result = dist (generator);
# else
    double result = gsl_ran_lognormal (rng.gsl_generator, meanlog, stdlog);
# endif
//     util::streamValidate(result);
    return result;
}

double random::max_multi_log_normal (double start, int n, double meanlog, double stdlog){
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
    double normp = pow( random::uniform_01(), 1.0 / n );
    double zval = gsl_cdf_ugaussian_Pinv (normp);
    double multi_sample = exp(meanlog + stdlog * zval);
    result = std::max(result, multi_sample);
# endif
//     util::streamValidate(result);
    return result;
}

double random::beta (double a, double b){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::beta_distribution<> dist (a, b);
    double result = dist(generator);
# else
    double result = gsl_ran_beta (rng.gsl_generator,a,b);
# endif
//     util::streamValidate(result);
    return result;
}
double random::betaWithMean (double m, double b){
    //TODO(performance): could do this calculation externally, and feed in a,b instead of mean,b
    double a = m * b / (1.0 - m);
//     util::streamValidate(a);
    return beta(a,b);
}

int random::poisson(double lambda){
    if( !(boost::math::isfinite)(lambda) ){
	//This would lead to an inifinite loop
	throw TRACED_EXCEPTION( "lambda is inf", Error::InfLambda );
    }
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::poisson_distribution<> dist (lambda);
    int result = dist(generator);
# else
    int result = gsl_ran_poisson (rng.gsl_generator, lambda);
# endif
//     util::streamValidate(result);
    return result;
}

bool random::bernoulli(double prob){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::bernoulli_distribution<> dist (prob);
    bool result = dist(generator);
# else
    assert( (boost::math::isfinite)(prob) );
    // return true iff our variate is less than the probability
    bool result =random::uniform_01() < prob;
# endif
//     util::streamValidate(result);
    return result;
}

int random::uniform(int n){
    assert( (boost::math::isfinite)(n) );
    return static_cast<int>( random::uniform_01() * n );
}

double random::weibull(double lambda, double k){
# ifdef OM_RANDOM_USE_BOOST_DIST
    boost::random::weibull_distribution<> dist (k, lambda);
    return dist(generator);
# else
    return gsl_ran_weibull( rng.gsl_generator, lambda, k );
# endif
}

} }
