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

#ifndef OM_util_random
#define OM_util_random

// Define to use boost distributions. Unfortunately these are not necessarily
// value-stable or portable.
// #define OM_RANDOM_USE_BOOST_DIST

#include "Global.h"
#include "util/errors.h"
#include <set>
#include <pcg_random.hpp>
#include <boost/random/uniform_01.hpp>
#include <gsl/gsl_rng.h>

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

namespace OM { namespace util {

// Support functions
template<class T>
long unsigned int sample_ulong (void *ptr) {
    T *rng = reinterpret_cast<T*>(ptr);
    BOOST_STATIC_ASSERT (sizeof(uint32_t) <= sizeof(long unsigned int));
    return static_cast<long unsigned int> ((*rng)());
}
template<class T>
double sample_double01 (void *ptr) {
    T *rng = reinterpret_cast<T*>(ptr);
    boost::random::uniform_01<T&> dist (*rng);
    return dist ();
}

template<class T>
gsl_rng_type make_gsl_rng_type(T& rng) {
    return {
        "OM_RNG",		// name
        rng.max(),		// max value
        rng.min(),		// min value
        0,				// size of state; not used here
        nullptr,			// re-seed function; don't use
        &sample_ulong<T>,
        &sample_double01<T>
    };
}

/// Our random number generator.
template<class T>
struct RNG {
    ///@brief Construction and checkpointing
    //@{
    /// Default constructor
    RNG(): m_rng(0) {
        m_gsl_type = make_gsl_rng_type(m_rng);
        m_gsl_gen.type = &m_gsl_type;
        m_gsl_gen.state = reinterpret_cast<void*>(&m_rng);
    }
    
    // Disable copying
    RNG(const RNG&) = delete;
    RNG& operator=(const RNG&) = delete;
  
    /// Allow moving, with explicit functions
    RNG(RNG&& other): m_rng(other.m_rng) {
        m_gsl_type = other.m_gsl_type;
        // The following point into instance of this RNG:
        m_gsl_gen.type = &m_gsl_type;
        m_gsl_gen.state = reinterpret_cast<void*>(&m_rng);
    }
    void operator=(RNG&& other) {
        m_rng = other.m_rng;
        m_gsl_type = other.m_gsl_type;
        // The following point into instance of this RNG:
        m_gsl_gen.type = &m_gsl_type;
        m_gsl_gen.state = reinterpret_cast<void*>(&m_rng);
    }
    
    /// Seed
    void seed(uint64_t seed) {
        m_rng.seed(seed);
    }
    
    /// Checkpointing
    /// 
    /// Note: this relies on the RNG having been constructed with the same seed.
    template<class S>
    void checkpoint(S& stream) {
        m_rng.binary_checkpoint(stream);
    }
    //@}
    
    ///@brief Random number distributions
    //@{
    /** Generate a random number in the range [0,1). */
    inline double uniform_01 () {
        boost::random::uniform_01<T&> rng_uniform01 (m_rng);
        return rng_uniform01 ();
    }
    
    /** This function returns a Gaussian random variate, with mean mean and
     * standard deviation std. The sampled value x ~ N(mean, std^2) . */
    double gauss (double mean, double std){
# ifdef OM_RANDOM_USE_BOOST_DIST
        boost::random::normal_distribution<> dist (mean, std);
        return dist(m_rng);
# else
        return gsl_ran_gaussian(&m_gsl_gen,std)+mean;
# endif
    }
    
    /** This function returns a random variate from the gamma distribution. */
    double gamma (double a, double b){
# ifdef OM_RANDOM_USE_BOOST_DIST
        boost::random::gamma_distribution<> dist (a, b);
        return dist(m_rng);
# else
        return gsl_ran_gamma(&m_gsl_gen, a, b);
# endif
    }
    
    /** This function returns a random variate from the lognormal distribution.
     * 
     * Mean is log(mean) - (sigma^2)/2.
     * Variance is (exp(sigma^2) - 1)*exp(2mu+sigma^2).
     * 
     * @param mu mean-log
     * @param sigma sigma-log
     */
    double log_normal (double meanlog, double stdlog){
# ifdef OM_RANDOM_USE_BOOST_DIST
        boost::random::lognormal_distribution<> dist (meanlog, stdlog);
        return dist (m_rng);
# else
        return gsl_ran_lognormal (&m_gsl_gen, meanlog, stdlog);
# endif
    }
    
    /** Return the maximum over multiple log-normal samples.
     *
     * @param start Initial value for running maximum (returns max of this and
     *      all samples)
     * @param n Number of samples
     * @param meanlog mean of underlying Gaussian
     * @param stdlog standard deviation of underlying Gaussian
     */
    double max_multi_log_normal (double start, int n, double meanlog, double stdlog){
# ifdef OM_RANDOM_USE_BOOST_DIST
        // We don't have a CDF available, so take log-normal samples
        double result = start;
        for (int i = 0; i < n; ++i) {
            double sample = log_normal(meanlog, stdlog);
            if (sample > result) result = sample;
        }
        return result;
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
        return std::max(start, multi_sample);
# endif
    }
    
    /** This function returns a random variate from the beta distribution. */
    double beta(double a, double b){
# ifdef OM_RANDOM_USE_BOOST_DIST
        boost::random::beta_distribution<> dist (a, b);
        return dist(m_rng);
# else
        return gsl_ran_beta (&m_gsl_gen,a,b);
# endif
    }
    
    /** This function wraps beta(), setting b=b and a such that m is the mean
     * of the distribution. */
    inline double betaWithMean(double m, double b) {
        //TODO(performance): could do this calculation externally, and feed in a,b instead of mean,b
        double a = m * b / (1.0 - m);
        return beta(a,b);
    }
    
    /** This function returns a random integer from the Poisson distribution with mean lambda. */
    int poisson(double lambda){
        if( !(boost::math::isfinite)(lambda) ){
            //This would lead to an inifinite loop
            throw TRACED_EXCEPTION( "lambda is inf", Error::InfLambda );
        }
# ifdef OM_RANDOM_USE_BOOST_DIST
        boost::random::poisson_distribution<> dist (lambda);
        return dist(m_rng);
# else
        return gsl_ran_poisson (&m_gsl_gen, lambda);
# endif
    }

    /** This function returns true with probability prob or 0 with probability
     * 1-prob (Bernoulli distribution). */
    bool bernoulli(double prob){
# ifdef OM_RANDOM_USE_BOOST_DIST
        boost::random::bernoulli_distribution<> dist (prob);
        return dist(m_rng);
# else
        assert( (boost::math::isfinite)(prob) );
        // return true iff our variate is less than the probability
        return uniform_01() < prob;
# endif
    }
    
    /** This function returns an integer from 0 to 1-n, where every value has
     * equal probability of being sampled. */
    inline int uniform (int n) {
        assert( (boost::math::isfinite)(n) );
        return static_cast<int>( uniform_01() * n );
    }
    
    /**
     * Return a variate sampled from the Weibull distribution.
     * 
     * The PDF is k · x^{k-1} exp{-(x/λ)^k} / λ^k
     * 
     * @param lambda (λ) is the scale parameter
     * @param k is the shape parameter
     */
    double weibull( double lambda, double k ){
# ifdef OM_RANDOM_USE_BOOST_DIST
        boost::random::weibull_distribution<> dist (k, lambda);
        return dist(m_rng);
# else
        return gsl_ran_weibull( &m_gsl_gen, lambda, k );
# endif
    }
    //@}
    
private:
    T m_rng;
    
    // Hooks for GSL distributions
    gsl_rng_type m_gsl_type;
    gsl_rng m_gsl_gen;
};

/// The global RNG
// I would prefer to use pcg64, but MSVC mysteriously fails
extern RNG<pcg32> global_RNG;

} }
#endif
