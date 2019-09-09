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

#include "Global.h"
#include <set>
#include <pcg_random.hpp>
#include <boost/random/uniform_01.hpp>
#include <gsl/gsl_rng.h>

namespace OM { namespace util {

/// Our random number generator.
struct RNG {
    ///@brief Construction and checkpointing
    //@{
    /// Default constructor
    RNG();
    /// Construction, given a seed
    explicit RNG(uint32_t seed);
    /// Construction from checkpoint
    explicit RNG(istream& stream);
    
    // Disable copying
    RNG(const RNG&) = delete;
    RNG& operator=(const RNG&) = delete;
  
    /// Allow moving, with explicit functions
    RNG(RNG&& other);
    void operator=(RNG&& other);
    
    /// Checkpointing
    void checkpoint (ostream& stream);
    //@}
    
    ///@brief Random number distributions
    //@{
    /** Generate a random number in the range [0,1). */
    inline double uniform_01 () {
        boost::random::uniform_01<pcg32&> rng_uniform01 (m_rng);
        return rng_uniform01 ();
    }
    
    /** This function returns a Gaussian random variate, with mean mean and
     * standard deviation std. The sampled value x ~ N(mean, std^2) . */
    double gauss (double mean, double std);
    
    /** This function returns a random variate from the gamma distribution. */
    double gamma (double a, double b);
    
    /** This function returns a random variate from the lognormal distribution.
     * 
     * Mean is log(mean) - (sigma^2)/2.
     * Variance is (exp(sigma^2) - 1)*exp(2mu+sigma^2).
     * 
     * @param mu mean-log
     * @param sigma sigma-log
     */
    double log_normal (double meanlog, double stdlog);
    
    /** Return the maximum over multiple log-normal samples.
     *
     * @param start Initial value for running maximum (returns max of this and
     *      all samples)
     * @param n Number of samples
     * @param meanlog mean of underlying Gaussian
     * @param stdlog standard deviation of underlying Gaussian
     */
    double max_multi_log_normal (double start, int n, double meanlog, double stdlog);
    
    /** This function returns a random variate from the beta distribution. */
    double beta(double a, double b);
    
    /** This function wraps beta(), setting b=b and a such that m is the mean
     * of the distribution. */
    inline double betaWithMean(double m, double b) {
        //TODO(performance): could do this calculation externally, and feed in a,b instead of mean,b
        double a = m * b / (1.0 - m);
        return beta(a,b);
    }
    
    /** This function returns a random integer from the Poisson distribution with mean lambda. */
    int poisson(double lambda);

    /** This function returns true with probability prob or 0 with probability
     * 1-prob (Bernoulli distribution). */
    bool bernoulli(double prob);
    
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
    double weibull( double lambda, double k );
    //@}
    
private:
    // I would prefer to use pcg64, but MSVC mysteriously fails
    pcg32 m_rng;
    
    // Hooks for GSL distributions
    gsl_rng_type m_gsl_type;
    gsl_rng m_gsl_gen;
};

/// The global RNG
extern RNG global_RNG;

} }
#endif
