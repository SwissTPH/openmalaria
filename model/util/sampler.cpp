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

#include "util/sampler.h"
#include "util/errors.h"
#include "util/random.h"
#include <cmath>

namespace OM { namespace util {

double NormalSample::asNormal( double mu, double sigma )const{
    return sigma*x + mu;
}

double NormalSample::asLognormal( double mu, double sigma )const{
    return exp( sigma*x + mu );
}

NormalSample NormalSample::generate(LocalRng& rng) {
    return NormalSample( rng.gauss(0.0, 1.0) );
}

NormalSample NormalSample::generate_correlated(NormalSample base, double correlation, double factor, LocalRng& rng) {
    if( correlation == 1.0 ) { return base; }
    
    double e = rng.gauss(0.0, factor);
    return NormalSample( base.x * correlation + e );
}

void NormalSampler::setParams( double m, double s ){
    mu = m;
    sigma = s;
}

void NormalSampler::setParams(const scnXml::SampledValueN& elt){
    mu = elt.getMean();
    if( elt.getDistr() == "const" ){
        if( elt.getSD().present() && elt.getSD().get() != 0.0) {
            throw util::xml_scenario_error( "attribute SD must be zero or omitted when distr=\"const\" or is omitted" );
        }
        sigma = 0.0;
        return;
    }
    
    if( elt.getDistr() != "normal" ){
        throw util::xml_scenario_error( "expected distr to be one of \"const\", \"lognormal\" (note: not all distributions are supported here)" );
    }
    if( !elt.getSD().present() ){
        throw util::xml_scenario_error( "attribute \"SD\" required for sampled normal value when distr is not const" );
    }
    sigma = elt.getSD().get();
}

double NormalSampler::sample(LocalRng& rng) const{
    if( sigma == 0.0 ){
        return mu;
    }
    return rng.gauss( mu, sigma );
}

void LognormalSampler::setParams( const scnXml::SampledValueLN& elt ){
    const double mean = elt.getMean();
    setParams( mean, elt );
}

void LognormalSampler::setParams( double mean, const scnXml::SampledValueCV& elt ){
    if( elt.getDistr() == "const" ){
        if( elt.getCV().present() && elt.getCV().get() != 0.0 ){
            throw util::xml_scenario_error( "lognormal distr: attribute CV must be zero or omitted when distr=\"const\" or is omitted" );
        }
        if( mean == 0.0 ){
            mu = -numeric_limits<double>::infinity();
        } else {
            mu = log(mean);
        }
        sigma = 0.0;
        return;
    }
    if( !elt.getCV().present() && !elt.getVariance().present())
        throw util::xml_scenario_error( "lognormal distr: attribute \"CV\" or \"variance\" required for sampled value when distr is not const" );
    if( elt.getCV().present() && elt.getVariance().present())
        throw util::xml_scenario_error( "lognormal distr: only one attribute \"CV\" or \"variance\" can be used for sampled value when distr is not const" );
    if( elt.getDistr() == "lognormal" ){
        if(elt.getCV().present())
            setMeanCV(mean, elt.getCV().get());
        else
            setMeanVariance(mean, elt.getVariance().get());
    }else{
        throw util::xml_scenario_error( "lognormal distr: expected distr to be one of \"const\", \"gamma\" (note: not all distributions are supported here)" );
    }
}

void LognormalSampler::setMeanCV( double mean, double CV ){
    this->CV = CV;

    // The distribution is "const"
    if( CV == 0.0 ){
        sigma = 0.0;
        // as a special case, we can support mean == CV == 0
        if( mean == 0.0 )
            mu = -numeric_limits<double>::infinity();
        else if(mean > 0)
            mu = log(mean);
        else
            throw util::xml_scenario_error( "const: required mean >= 0" );
        return;
    }

    // If the distirbution is not "const"
    if( mean < 0 )
        throw util::xml_scenario_error( "log-normal: required mean > 0" );
    if( CV < 0 )
        throw util::xml_scenario_error( "log-normal: required CV >= 0" );
    
    // The following formulae are derived from:
    // X ~ lognormal(μ, σ)
    // E(X) = exp(μ + σ² / 2)
    // Var(X) = exp(σ² + 2μ)(exp(σ²) - 1)
    // Var = (CV * mean)²
    
    const double a = 1 + (CV * CV);
    mu = log( mean / sqrt(a) );
    sigma = sqrt(log(a));
}

void LognormalSampler::setMeanVariance( double mean, double variance ){
    if( mean <= 0 )
        throw util::xml_scenario_error( "log-normal: required mean > 0" );
    if( variance < 0 )
        throw util::xml_scenario_error( "log-normal: required variance >= 0" );

    // The following formulae are derived from:
    // X ~ lognormal(μ, σ)
    // E(X) = exp(μ + σ² / 2)
    // Var(X) = exp(σ² + 2μ)(exp(σ²) - 1)
    // Var = (CV * mean)²

    if( variance == 0.0)
    {
        mu = log(mean);
        sigma = 0.0;
        return;
    }

    const double CV = variance / mean;
    const double a = 1 + (CV * CV);
    mu = log( mean / sqrt(a) );
    sigma = sqrt(log(a));
}

void LognormalSampler::scaleMean(double scalar){
    mu += log(scalar);
}

double LognormalSampler::mean() const{
    return exp(mu + 0.5*sigma*sigma);
}

double LognormalSampler::sample(LocalRng& rng) const{
    if( sigma == 0.0 ){
        return exp( mu );
    } else {
        return rng.log_normal( mu, sigma );
    }
}

double LognormalSampler::cdf(double x) const {
    // 1) anything at or below 0 has CDF = 0
    if (x <= 0.0)
        return 0.0;

    // 2) Degenerate case
    if (sigma == 0.0)
    {
        if(log(x) >= mu) return 1.0;
        else return 0.0;
    }
    return gsl_cdf_lognormal_P(x, mu, sigma);
}

void GammaSampler::setParams( double mean, const scnXml::SampledValueCV& elt ){
    if( elt.getDistr() == "const" ){
        if( elt.getCV().present() && elt.getCV().get() != 0.0 ){
            throw util::xml_scenario_error( "gamma distr: attribute CV must be zero or omitted when distr=\"const\" or is omitted" );
        }
        if( mean == 0.0 ){
            mu = -numeric_limits<double>::infinity();
        } else {
            mu = mean;
        }
        return;
    }
    if( !elt.getCV().present() && !elt.getVariance().present())
        throw util::xml_scenario_error( "gamma distr: attribute \"CV\" or \"variance\" required for sampled value when distr is not const" );
    if( elt.getCV().present() && elt.getVariance().present())
        throw util::xml_scenario_error( "gamma distr: only one attribute \"CV\" or \"variance\" can be used for sampled value when distr is not const" );
    if( elt.getDistr() == "gamma" ){
        if(elt.getCV().present())
            setMeanCV(mean, elt.getCV().get());
        else
            setMeanVariance(mean, elt.getVariance().get());
    }else{
        throw util::xml_scenario_error( "gamma distr: expected distr to be one of \"const\", \"gamma\" (note: not all distributions are supported here)" );
    }
}

void GammaSampler::setMeanCV( double mean, double CV ){
    if( mean <= 0 )
        throw util::xml_scenario_error( "gamma: required mean > 0" );
    if( CV < 0 )
        throw util::xml_scenario_error( "gamma: required CV >= 0" );

    mu = mean;
    this->CV = CV;

    if( CV == 0.0 )
        return;

    // 1 / sqrt(k) = CV
    // sqrt(k) = 1/CV
    // k = 1 / CV^2
    k = 1.0/(CV*CV);
    // k * theta = mean
    // theta = mean / k
    theta = mu / k;
}

void GammaSampler::setMeanVariance( double mean, double variance ){
    mu = mean;

    if( mean <= 0 )
        throw util::xml_scenario_error( "gamma: required mean > 0" );
    if( variance < 0 )
        throw util::xml_scenario_error( "gamma: required variance >= 0" );

    if( variance == 0.0 )
        return;

    this->variance = variance;
    // sigma / mu = 1 / sqrt(k)
    // sqrt(k) = mu / sigma
    // k = mu^2 / variance
    k = (mu*mu)/this->variance;
    // k * theta = mean
    // theta = mean / k
    theta = mu / k;
}

void GammaSampler::scaleMean(double scalar) {
    if (scalar <= 0.0) {
        // Invalid scalar, return without making changes
        return;
    }

    // Scale the mean
    mu *= scalar;

    if (!std::isnan(this->CV) && this->CV > 0.0) {
        this->variance = (mu * this->CV) * (mu * this->CV);
    }

    if (!std::isnan(this->variance) && this->variance > 0.0) {
        // Recalculate k and theta based on the scaled mean and fixed variance
        k = (mu * mu) / this->variance;  // Shape parameter
        theta = this->variance / mu;    // Scale parameter
    }
}

double GammaSampler::mean() const {
    return mu;
}

double GammaSampler::sample(LocalRng& rng) const{
    if(isnan(theta)) // CV=0
        return mu;
    return rng.gamma(k, theta);
}

double GammaSampler::cdf(double x) const {
    // 1) anything at or below 0 has CDF = 0
    if (x <= 0.0)
        return 0.0;

    // 2) Degenerate case
    if(isnan(theta)) // CV=0 or variance=0
    {
        if(x >= mu) return 1.0;
        else return 0.0;
    }

    return gsl_cdf_gamma_P(x, k, theta);
}
        
void BetaSampler::setParamsMV( double mean, double variance ){
    if( variance > 0.0 ){
        // double c = mean / (1.0 - mean);
        // double cp1 = c + 1.0;
        // double s = cp1*cp1*variance;
        // b = (s - c) / (s*cp1);
        // a = c*b;
        
        if (variance >= mean * (1.0 - mean))
            throw util::xml_scenario_error("BetaSampler::setParamsMV: require variance < mean(1.0 - mean)");

        double c = mean * ((1.0 - mean) / variance) - 1.0;
        a = mean * c;
        b = c - a;
    }else if (variance == 0.0)
    {
        // Not using the Beta distirbution, we simply return the mean when sampling
        if( variance != 0.0 ){
            throw util::xml_scenario_error("BetaSampler::setParamsMV: require variance ≥ 0");
        }
        a = mean;
        b = 0.0;
    }
    else
        throw util::xml_scenario_error("BetaSampler::setParamsMV: require variance ≥ 0");
}

double BetaSampler::sample(LocalRng& rng) const{
    if( b == 0.0 ){
        return a;
    }else{
        assert( a>0.0 && b>0.0 );
        return rng.beta( a, b );
    }
}


} }
