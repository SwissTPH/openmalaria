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

#include "Global.h"
#include <set>

namespace OM { namespace util {

/** Random number generator.
 *
 * This interface should be independant of implementation. */
namespace random {
    ///@brief Setup & cleanup; checkpointing
    //@{
    /// Reseed the random-number-generator with seed (usually InputData.getISeed()).
    void seed (uint32_t seed);
    
    void checkpoint (istream& stream, int seedFileNumber);
    void checkpoint (ostream& stream, int seedFileNumber);
    //@}
    
    ///@brief Random number distributions
    //@{
    /** Generate a random number in the range [0,1). */
    double uniform_01 ();
    
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
    double log_normal (double mu, double sigma);
    
    /** Used for performance reasons. Calling rngLogNormal 5 times is 50% slower. */
    double sampleFromLogNormal (double normp, double meanlog, double stdlog);
    
    /** This function returns a random variate from the beta distribution. */
    double beta(double a, double b);
    
    /** This function wraps beta(), setting b=b and a such that m is the mean
     * of the distribution. */
    double betaWithMean(double m, double b);
    
    /** This function returns a random integer from the Poisson distribution with mean lambda. */
    int poisson(double lambda);

    /** This function returns true with probability prob or 0 with probability
     * 1-prob (Bernoulli distribution). */
    bool bernoulli(double prob);
    
    /** This function returns an integer from 0 to 1-n, where every value has
     * equal probability of being sampled. */
    int uniform (int n);
    
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
}
} }
