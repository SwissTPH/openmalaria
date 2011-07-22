/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include "schema/util.h"
#include <limits>

namespace OM { namespace util {
    
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
        
        static NormalSample generate();
        
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
        inline void setParams( const scnXml::NormalSample& elt ){
            setParams( elt.getMu(), elt.getSigma() );
        }
        
        /** Sample a value. */
        double sample() const;
        
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
    class LognormalSampler {
    public:
        LognormalSampler() :
            mu( numeric_limits<double>::signaling_NaN() ),
            sigma( numeric_limits<double>::signaling_NaN() )
        {}
        
        /** Set parameters such that samples taken are:
         * X ~ log N( log(mean)-s²/2, s² )
         * 
         * @param mean Mean of sampled variates
         * @param s Square-root of variance of logarithm of sampled variates
         */
        void setParams( double mean, double s );
        /** Set the mean, leave sigma unchanged. */
        void setMean( double mean );
        /** As above, using an XML element. */
        inline void setParams( const scnXml::LognormalSample& elt ){
            setParams( elt.getMean(), elt.getSigma() );
        }
        
        /** Sample a value. */
        double sample() const;
        
        /** Create a log-normal sample from an existing normal sample. */
        inline double sample(NormalSample sample) const{
            return sample.asLognormal( mu, sigma );
        }
        
    private:
        double mu, sigma;
    };
    
    /** Sampler for beta distribution.
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
        double sample() const;
        
    private:
        //Note: if b is 0, then alpha is mean. Otherwise a and b are
        //the expected (α,β) parameters to the beta distribution.
        double a, b;
    };
    
} }
#endif

