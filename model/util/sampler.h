/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
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

#ifndef UTIL_SAMPLER
#define UTIL_SAMPLER

#include "Global.h"
#include "util/random.h"
#include "schema/util.h"
#include <limits>

namespace OM { namespace util {
    
    class Sampler
    {
    public:
        virtual ~Sampler() = default;
        virtual double sample(LocalRng& rng) const = 0;
        virtual void scaleMean( double scalar ) = 0;
        virtual double mean() const = 0;
        virtual double cdf(double x) const { throw std::runtime_error("cdf() not implemented for this distribution"); }
    };

    /** A normal sample, which can be turned into various log-normal samples.
     * 
     * Allows generation of correlated log-normal samples with different sigma.
     */
    class NormalSample {
    public:
        NormalSample() : x(numeric_limits<double>::signaling_NaN()) {}
        
        /// convert sample to N(mu,sigma)
        double asNormal( double mu, double sigma ) const;

        /// convert sample to lnN(mu,sigma)
        double asLognormal( double mu, double sigma ) const;
        
        static NormalSample generate(LocalRng& rng);
        
        /// Generate a sample correlated with base.
        /// 
        /// @param base sample correlated against
        /// @param correlation expected correlation between these parameters
        ///     (note: if this is used to sample a log-normal, then the
        ///     correlation is on the log-scale)
        /// @param factor should equal sqrt(1 - correlation^2); may be cached
        static NormalSample generate_correlated(NormalSample base, double correlation, double factor, LocalRng& rng);
        
    private:
        NormalSample( double variate ) : x(variate) {}
        
        double x;	// variate (sampled from N(0,1))
    };
    
    /** Sampler for normal values */
    class NormalSampler {
    public:
        NormalSampler() :
            mu( numeric_limits<double>::signaling_NaN() ),
            sigma( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters such that samples taken are:
         * X ~ N( m, s² )
         * 
         * @param m Mean of sampled variates
         * @param s Square-root of variance of sampled variates
         */
        void setParams( double m, double s );

        /** As above, using an XML element. */
        void setParams( const scnXml::SampledValueN& elt );
        
        /** Sample a value. */
        double sample(LocalRng& rng) const;

        /** Create a log-normal sample from an existing normal sample. */
        inline double sample(NormalSample sample) const{
            return sample.asNormal( mu, sigma );
        }
        
        /// Return mu / mean of distribution
        inline double getMu() const{
            return mu;
        }

        /// Return sigma / standard deviation of distribution
        inline double getSigma() const{
            return sigma;
        }
        
    private:
        double mu, sigma;
    };
    
    /** Sampler for log-normal values */
    class LognormalSampler : public Sampler {
    public:
        LognormalSampler() :
            mu( numeric_limits<double>::signaling_NaN() ),
            sigma( numeric_limits<double>::signaling_NaN() ),
            CV( numeric_limits<double>::signaling_NaN() )
        {}

        LognormalSampler(double mean, const scnXml::SampledValueCV& elt) :
            mu( numeric_limits<double>::signaling_NaN() ),
            sigma( numeric_limits<double>::signaling_NaN() )
        {
            setParams(mean, elt);
        }
        
        /// Set specified mean and CV from XML element
        void setParams( const scnXml::SampledValueLN& elt );

        /// Set specified mean and CV from XML element
        void setParams( double mean, const scnXml::SampledValueCV& elt );

        /** Set log-normal parameters from mean and CV. */
        void setMeanCV( double mean, double CV );

        /** Set log-normal parameters from mean and variance. */
        void setMeanVariance( double mean, double CV );

        /** Scale the mean (i.e. multiply by a scalar).
         * 
         * Note that sigma (when specified via CV) is independent of the mean,
         * so one can safely multiply the mean without affecting CV. */
        void scaleMean( double scalar );
        
        /** Get the mean. */
        double mean() const;

        /** Sample a value. */
        double sample(LocalRng& rng) const;
        
        /**
         * Compute the cumulative distribution function (CDF) of this log-normal sampler.
         *
         * Special cases:
         *  - If x ≤ 0, returns 0.
         *  - If sigma == 0, the distribution collapses to a point mass,
         * so returns 1 if log(x) ≥ mu, else 0.
         *  - Otherwise uses gsl_cdf_lognormal_P(x, mu, sigma).
         *
         * @param x  The point at which to evaluate the CDF (must be a real number).
         * @return   P(X ≤ x), where X ~ LogNormal(mu, sigma²).
         */
        double cdf(double x) const;

        /** Create a log-normal sample from an existing normal sample. */
        inline double sample(NormalSample sample) const{
            return sample.asLognormal( mu, sigma );
        }
        
        /** Return true if and only if parameters have been set. */
        inline bool isSet() const{
            return mu == mu;    // mu is NaN iff not set
        }
        
    private:
        // log-space parameters
        double mu, sigma;
        double CV;
    };

    /** Sampler for gamma values */
    class GammaSampler : public Sampler {
    public:
        GammaSampler(const scnXml::SampledValueCV& elt) :
            mu( 1.0 ),
            k( numeric_limits<double>::signaling_NaN() ),
            theta( numeric_limits<double>::signaling_NaN() ),
            variance( numeric_limits<double>::signaling_NaN() ),
            CV( numeric_limits<double>::signaling_NaN() )
        {
            setParams(mu, elt);
        }
        
        /// Set specified mean and CV from XML element
        void setParams( double mean, const scnXml::SampledValueCV& elt );
        
        /** Set gamma parameters from mean and CV. */
        void setMeanCV( double mean, double CV );

        /** Set gamma parameters from mean and variance. */
        void setMeanVariance( double mean, double variance );

        /** Scale the mean (i.e. multiply by a scalar).
         * 
         * Note that sigma (when specified via CV) is independent of the mean,
         * so one can safely multiply the mean without affecting CV. */
        void scaleMean( double scalar );
        
        /** Get the mean. */
        double mean() const;

        /** Sample a value. */
        double sample(LocalRng& rng) const;
        
        /**
         * Compute the cumulative distribution function (CDF) of this gamma sampler.
         *
         * Special cases:
         *  - If x ≤ 0, returns 0.
         *  - If CV == 0 or variance == 0, the distribution collapses to a point mass,
         * so returns 1 if x ≥ mu, else 0.
         *  - Otherwise uses gsl_cdf_gamma_P(x, k, theta) where
         *      • k     = shape parameter
         *      • theta = scale parameter
         *
         * @param x  The point at which to evaluate the CDF (must be a real number).
         * @return   P(X ≤ x), where X ~ Gamma(k, θ).
         */
        double cdf(double x) const;
        
    private:
        double mu, k, theta;
        double variance;
        double CV;
    };
    
    /** Sampler for the Beta distribution.
     *
     * Input may be alpha and beta or mean and variance. As a special case,
     * variance zero is also supported, meaning return the mean. */
    class BetaSampler {
    public:
        BetaSampler() :
            a( numeric_limits<double>::signaling_NaN() ),
            b( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters: alpha, beta.
         * 
         * As a special case, when beta=0, alpha is interpreted as the mean and
         * is returned without sampling. */
        inline void setParams( double alpha, double beta ){
            a=alpha;
            b=beta;
        }
        
        /** Set parameters: alpha and beta are calculated such that the mean
         * and variance are as given. Variance may be zero, in which case the
         * mean is returned directly without sampling. */
        void setParamsMV( double mean, double variance );
        
        /** Set parameters from an XML element. */
        void setParams( const scnXml::BetaMeanSample& elt ){
            setParamsMV( elt.getMean(), elt.getVariance() );
        }
        
        /** Sample a value. */
        double sample(LocalRng& rng) const;
        
    private:
        //Note: if b is 0, then alpha is mean. Otherwise a and b are
        //the expected (α,β) parameters to the beta distribution.
        double a, b;
    };
    
    /** Sampler for the Weibull distribution. */
    class WeibullSampler {
    public:
        WeibullSampler() :
            scale( numeric_limits<double>::signaling_NaN() ),
            shape( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters from scale, shape. */
        inline void setScaleShape( double lambda, double k ){
            scale = lambda;
            shape = k;
        }
        
        /** Set parameters from an XML element. */
        inline void setParams( const scnXml::WeibullSample& elt ){
            setScaleShape( elt.getScale(), elt.getShape() );
        }
        
        /** Sample a value. */
        inline double sample(LocalRng& rng) const{
            return rng.weibull( scale, shape );
        }
        
    private:
        double scale, shape;    // λ, k
    };
    
} }
#endif
