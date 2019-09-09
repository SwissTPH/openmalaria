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

NormalSample NormalSample::generate() {
    return NormalSample( global_RNG.gauss(0.0, 1.0) );
}

NormalSample NormalSample::generate_correlated(NormalSample base, double correlation, double factor) {
    if( correlation == 1.0 ) { return base; }
    
    double e = global_RNG.gauss(0.0, factor);
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
double NormalSampler::sample() const{
    if( sigma == 0.0 ){
        return mu;
    }
    return global_RNG.gauss( mu, sigma );
}

void LognormalSampler::setParams( const scnXml::SampledValueLN& elt ){
    const double mean = elt.getMean();
    setParams( mean, elt );
}

void LognormalSampler::setParams( double mean, const scnXml::SampledValueCV& elt ){
    if( elt.getDistr() == "const" ){
        if( elt.getCV().present() && elt.getCV().get() != 0.0 ){
            throw util::xml_scenario_error( "attribute CV must be zero or omitted when distr=\"const\" or is omitted" );
        }
        mu = log(mean);
        sigma = 0.0;
        return;
    }
    
    if( !elt.getCV().present() ){
        throw util::xml_scenario_error( "attribute \"CV\" required for sampled value when distr is not const" );
    }
    
    if( elt.getDistr() == "lognormal" ){
        setMeanCV( mean, elt.getCV().get() );
    }else{
        throw util::xml_scenario_error( "expected distr to be one of \"const\", \"lognormal\" (note: not all distributions are supported here)" );
    }
}

void LognormalSampler::setMeanCV( double mean, double CV ){
    if( CV == 0.0 ){
        sigma = 0.0;
        // as a special case, we can support mean == CV == 0
        if( mean == 0.0 ){
            mu = -numeric_limits<double>::infinity();
        } else {
            mu = log(mean);
        }
        return;
    }
    if( !(mean > 0.0 && CV > 0.0) ){
        throw util::xml_scenario_error( "log-normal: required mean > 0 and CV ≥ 0" );
    }
    
    // The following formulae are derived from:
    // X ~ lognormal(μ, σ)
    // E(X) = exp(μ + σ² / 2)
    // Var(X) = exp(σ² + 2μ)(exp(σ²) - 1)
    // Var = CV * mean
    
    const double a = 1 + CV * CV;
    mu = log( mean / sqrt(a) );
    sigma = sqrt(log(a));
}
void LognormalSampler::scaleMean(double scalar){
    mu += log(scalar);
}
double LognormalSampler::mean() const{
    return exp(mu + 0.5*sigma*sigma);
}
double LognormalSampler::sample() const{
    if( sigma == 0.0 ){
        return exp( mu );
    } else {
        return global_RNG.log_normal( mu, sigma );
    }
}

void BetaSampler::setParamsMV( double mean, double variance ){
    if( variance > 0.0 ){
        double c = mean / (1.0 - mean);
        double cp1 = c + 1.0;
        double s = cp1*cp1*variance;
        b = (s - c) / (s*cp1);
        a = c*b;
    }else{
        if( variance != 0.0 ){
            throw TRACED_EXCEPTION_DEFAULT("BetaSampler::setParamsMV: require variance ≥ 0");
        }
        a = mean;
        b = 0.0;
    }
}
double BetaSampler::sample() const{
    if( b == 0.0 ){
        return a;
    }else{
        assert( a>0.0 && b>0.0 );
        return global_RNG.beta( a, b );
    }
}


} }
