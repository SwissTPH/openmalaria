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

#ifndef Hmod_DecayFunction
#define Hmod_DecayFunction

#include <limits>
#include <memory>

#include "Global.h"
#include "util/sampler.h"
#include "util/UnitParse.h"

namespace OM {
namespace util {

/** An interface for a few types of decay function (some of which may also be
 * suitible survival functions).
 *
 * Heterogeneity is implemented by passing a DecayFunction object to the eval
 * function; this should be sampled by the same DecayFunction.
 *****************************************************************************/
template <typename SamplerT = Sampler>
class DecayFunction
{
public:
    DecayFunction(bool increasing = false, double initialEfficacy = 1.0, double CV = 0.0) : 
        increasing(increasing), initialEfficacy(initialEfficacy), het(std::move(LognormalSampler::fromMeanCV(1.0, CV))) {}

    virtual ~DecayFunction() {}
    
    /** Return a new decay function, constructed from an XML element.
     * 
     * @param elt XML element specifying which function to use and parameters
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    static unique_ptr<DecayFunction<SamplerT>> makeObject(
        const scnXml::DecayFunction& elt, const char* eltName
    );

    /** Sample a DecayFunction value (should be stored per individual).
     * 
     * Note that a DecayFunction is needed to call eval() even if heterogeneity
     * is not wanted. If sigma = 0 then the random number stream will not be
     * touched. */
    unique_ptr<DecayFunction<SamplerT>> hetSample(LocalRng& rng) const {
        return hetSample(het->sample(rng));
    }
    
    /** Generate a DecayFunction value from an existing sample. */
    template <typename SampleT>
    std::unique_ptr<DecayFunction<SamplerT>> hetSample(SampleT sample) const {
        return hetSample(het->sample(sample));
    }
    
    /** Generate a DecayFunction value from an existing sample. */
    virtual unique_ptr<DecayFunction<SamplerT>> hetSample(double hetFactor) const =0;

    /** Say you have a population of objects which each have two states:
     * decayed and not decayed. If you want to use a DecayFunction to model
     * the proportion of objects which have decayed, you need to work out per
     * object the age of decay. This function does that.
     * 
     * This is only valid where mu and sigma parameters are zero.
     * 
     * @returns Age at which an object should decay. */
    virtual SimTime sampleAgeOfDecay (LocalRng& rng) const =0;
    
    inline double eval(double ageDays) const {
        if(increasing)
            return 1.0 - compute( ageDays ) * initialEfficacy;
        else
            return compute( ageDays ) * initialEfficacy;
    }

    virtual double compute(double ageDays) const = 0;

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        // timeFactorHet & stream;
    }

private:
    bool increasing;
    double initialEfficacy;
    std::shared_ptr<SamplerT> het;
};

inline double readLToDays( const scnXml::DecayFunction& elt ){
    if( !elt.getL().present() ){
        throw xml_scenario_error( "decay function: attribute L required" );
    }
    return UnitParse::durationToDays(elt.getL().get(), UnitParse::YEARS);
}

template<typename SamplerT = Sampler>
class ConstantDecayFunction : public DecayFunction<SamplerT> {
public:
    ConstantDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction<SamplerT>(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const {
        // Note: we now require all decay functions to return 0 when time > 0
        // and the DecayFunction is default-constructed. So const *after deployment*.
        if( effectiveAge * hetFactor == numeric_limits<double>::infinity() )
            return 0.0;
        return 1.0;
    }
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return sim::future();        // decay occurs "in the future" (don't use sim::never() because that is interpreted as being in the past)
    }

    unique_ptr<DecayFunction<SamplerT>> hetSample(double hetFactor) const {
        unique_ptr<ConstantDecayFunction<SamplerT>> copy = std::make_unique<ConstantDecayFunction<SamplerT>>(*this);
        copy->hetFactor = hetFactor;
        return std::move(copy);
    }

private:
    double hetFactor;
};

template<typename SamplerT = Sampler>
class StepDecayFunction : public DecayFunction<SamplerT> {
public:
    StepDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction<SamplerT>(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        if( effectiveAge * invL * hetFactor < 1.0 ){
            return 1.0;
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( 1.0 / invL );
    }

    unique_ptr<DecayFunction<SamplerT>> hetSample(double hetFactor) const {
        unique_ptr<StepDecayFunction> copy = make_unique<StepDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, hetFactor;
};

template<typename SamplerT = Sampler>
class LinearDecayFunction : public DecayFunction<SamplerT> {
public:
    LinearDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction<SamplerT>(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        // cout << "Linear " << effectiveAge << " " << invL << " " << hetFactor <<endl;
        if( effectiveAge * invL * hetFactor < 1.0 ){
            return 1.0 - effectiveAge * invL * hetFactor;
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        // Note: rounds to nearest. Object may decay instantly or at time L.
        return sim::roundToTSFromDays(rng.uniform_01() / invL);
    }

    unique_ptr<DecayFunction<SamplerT>> hetSample(double hetFactor) const {
        unique_ptr<LinearDecayFunction> copy = make_unique<LinearDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, hetFactor;
};

template<typename SamplerT = Sampler>
class ExponentialDecayFunction : public DecayFunction<SamplerT> {
public:
    ExponentialDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction<SamplerT>(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( log(2.0) / readLToDays(elt) ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        return exp( -effectiveAge * invL * hetFactor);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( -log(rng.uniform_01()) / invL );
    }

    unique_ptr<DecayFunction<SamplerT>> hetSample(double hetFactor) const {
        unique_ptr<ExponentialDecayFunction> copy = make_unique<ExponentialDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, hetFactor;
};

template<typename SamplerT = Sampler>
class WeibullDecayFunction : public DecayFunction<SamplerT> {
public:
    WeibullDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction<SamplerT>(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        constOverLambda( pow(log(2.0),1.0/elt.getK()) / readLToDays(elt) ),
        k( elt.getK() ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        double p = -pow(effectiveAge * constOverLambda * hetFactor, k);
        if(p < -700.0)
            return 0.0;
        else
            return exp(p);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( pow( -log(rng.uniform_01()), 1.0/k ) / constOverLambda );
    }

    unique_ptr<DecayFunction<SamplerT>> hetSample(double hetFactor) const {
        unique_ptr<WeibullDecayFunction> copy = make_unique<WeibullDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double constOverLambda, k, hetFactor;
};

template<typename SamplerT = Sampler>
class HillDecayFunction : public DecayFunction<SamplerT> {
public:
    HillDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction<SamplerT>(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        k( elt.getK() ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        return 1.0 / (1.0 + pow(effectiveAge * invL * hetFactor, k));
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( pow( 1.0 / rng.uniform_01() - 1.0, 1.0/k ) / invL );
    }

    unique_ptr<DecayFunction<SamplerT>> hetSample(double hetFactor) const {
        unique_ptr<HillDecayFunction> copy = make_unique<HillDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, k, hetFactor;
};

template<typename SamplerT = Sampler>
class SmoothCompactDecayFunction : public DecayFunction<SamplerT> {
public:
    SmoothCompactDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction<SamplerT>(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( 1.0 / readLToDays(elt) ),
        k( elt.getK() ),
        hetFactor(0.0)
    {}

    double compute(double effectiveAge) const{
        if( effectiveAge * invL * hetFactor < 1.0 ){
            return exp( k - k / (1.0 - pow(effectiveAge * invL * hetFactor, 2.0)) );
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( sqrt( 1.0 - k / (k - log( rng.uniform_01() )) ) / invL );
    }

    unique_ptr<DecayFunction<SamplerT>> hetSample(double hetFactor) const {
        unique_ptr<SmoothCompactDecayFunction> copy = make_unique<SmoothCompactDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, k, hetFactor;
};

template<class T, typename SamplerT = Sampler>
class OperatorDecayFunction : public DecayFunction<SamplerT> {
public:
    OperatorDecayFunction( const scnXml::DecayFunction& elt ) : 
        DecayFunction<SamplerT>(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()) {
        const scnXml::DecayFunction::DecaySequence &decaySequence = elt.getDecay();
        if(decaySequence.size() != 2)
            throw xml_scenario_error("Operator decay function expects two decay functions, " + to_string(decaySequence.size()) +"  were given.");

        f1 = DecayFunction<SamplerT>::makeObject(decaySequence[0], "Operator::f1");
        f2 = DecayFunction<SamplerT>::makeObject(decaySequence[1], "Operator::f2");
    }

    OperatorDecayFunction(const OperatorDecayFunction &copy, unique_ptr<DecayFunction<SamplerT>> f1, unique_ptr<DecayFunction<SamplerT>> f2) : 
        DecayFunction<SamplerT>(copy),
        f1(move(f1)), f2(move(f2)) {}

    double compute(double effectiveAge) const {
        return max(min(op(f1->eval(effectiveAge), f2->eval(effectiveAge)), 1.0), 0.0);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return sim::roundToTSFromDays( max(f1->sampleAgeOfDecay(rng), f2->sampleAgeOfDecay(rng)) );
    }

    unique_ptr<DecayFunction<SamplerT>> hetSample(double hetFactor) const {
        unique_ptr<DecayFunction<SamplerT>> f1hetSample = f1->hetSample(hetFactor);
        unique_ptr<DecayFunction<SamplerT>> f2hetSample = f2->hetSample(hetFactor);
        unique_ptr<OperatorDecayFunction> copy = make_unique<OperatorDecayFunction>(*this, move(f1hetSample), move(f2hetSample));
        return move(copy);
    }

private:
    unique_ptr<DecayFunction<SamplerT>> f1, f2;
    T op;
};

// -----  interface / static functions  -----
template<typename SamplerT>
unique_ptr<DecayFunction<SamplerT>> DecayFunction<SamplerT>::makeObject(
    const scnXml::DecayFunction& elt, const char* eltName
){
    // Type mostly equivalent to a std::string:
    const scnXml::Function& func = elt.getFunction();
    if( func == "constant" ){
        return make_unique<ConstantDecayFunction<SamplerT>>( elt );
    }else if( func == "step" ){
        return make_unique<StepDecayFunction<SamplerT>>( elt );
    }else if( func == "linear" ){
        return make_unique<LinearDecayFunction<SamplerT>>( elt );
    }else if( func == "exponential" ){
        return make_unique<ExponentialDecayFunction<SamplerT>>( elt );
    }else if( func == "weibull" ){
        return make_unique<WeibullDecayFunction<SamplerT>>( elt );
    }else if( func == "hill" ){
        return make_unique<HillDecayFunction<SamplerT>>( elt );
    }else if( func == "smooth-compact" ){
        return make_unique<SmoothCompactDecayFunction<SamplerT>>( elt );
    }else if( func == "plus" ){
        return make_unique<OperatorDecayFunction<std::plus<double>, SamplerT>>( elt );
    }else if( func == "minus" ){
        return make_unique<OperatorDecayFunction<std::minus<double>, SamplerT>>( elt );
    }else if( func == "divides" ){
        return make_unique<OperatorDecayFunction<std::divides<double>, SamplerT>>( elt );
    }else if( func == "multiplies" ){
        return make_unique<OperatorDecayFunction<std::multiplies<double>, SamplerT>>( elt );
    }else{
        throw xml_scenario_error("decay function type " + string(func) + " of " + string(eltName) + " unrecognized");
    }
}


} }
#endif
