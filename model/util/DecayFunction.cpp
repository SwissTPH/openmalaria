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

#include "util/DecayFunction.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/random.h"
#include "util/UnitParse.h"

#include <cmath>
#include <functional>
#include <stdexcept>

namespace OM {
namespace util {

double readLToDays( const scnXml::DecayFunction& elt ){
    if( !elt.getL().present() ){
        throw util::xml_scenario_error( "decay function: attribute L required" );
    }
    return UnitParse::durationToDays(elt.getL().get(), UnitParse::YEARS);
}

class ConstantDecayFunction : public DecayFunction {
public:
    ConstantDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
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

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<ConstantDecayFunction> copy = make_unique<ConstantDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }

private:
    double hetFactor;
};

class StepDecayFunction : public DecayFunction {
public:
    StepDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
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

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<StepDecayFunction> copy = make_unique<StepDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, hetFactor;
};

class LinearDecayFunction : public DecayFunction {
public:
    LinearDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
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

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<LinearDecayFunction> copy = make_unique<LinearDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, hetFactor;
};

class ExponentialDecayFunction : public DecayFunction {
public:
    ExponentialDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
        invL( log(2.0) / readLToDays(elt) ),
        hetFactor(0.0)
    {}
    
    double compute(double effectiveAge) const{
        return exp( -effectiveAge * invL * hetFactor);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( -log(rng.uniform_01()) / invL );
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<ExponentialDecayFunction> copy = make_unique<ExponentialDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, hetFactor;
};

class WeibullDecayFunction : public DecayFunction {
public:
    WeibullDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
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

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<WeibullDecayFunction> copy = make_unique<WeibullDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double constOverLambda, k, hetFactor;
};

class HillDecayFunction : public DecayFunction {
public:
    HillDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
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

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<HillDecayFunction> copy = make_unique<HillDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, k, hetFactor;
};

class SmoothCompactDecayFunction : public DecayFunction {
public:
    SmoothCompactDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()),
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

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<SmoothCompactDecayFunction> copy = make_unique<SmoothCompactDecayFunction>(*this);
        copy->hetFactor = hetFactor;
        return move(copy);
    }
    
private:
    double invL, k, hetFactor;
};

template<class T>
class OperatorDecayFunction : public DecayFunction {
public:
    OperatorDecayFunction( const scnXml::DecayFunction& elt ) : 
        DecayFunction(elt.getIncreasing(), elt.getInitialEfficacy(), elt.getCV()) {
        const scnXml::DecayFunction::DecaySequence &decaySequence = elt.getDecay();
        if(decaySequence.size() != 2)
            throw util::xml_scenario_error("Operator decay function expects two decay functions, " + to_string(decaySequence.size()) +"  were given.");

        f1 = makeObject(decaySequence[0], "Operator::f1");
        f2 = makeObject(decaySequence[1], "Operator::f2");
    }

    OperatorDecayFunction(const OperatorDecayFunction &copy, unique_ptr<DecayFunction> f1, unique_ptr<DecayFunction> f2) : 
        DecayFunction(copy),
        f1(move(f1)), f2(move(f2)) {}

    double compute(double effectiveAge) const {
        return max(min(op(f1->eval(effectiveAge), f2->eval(effectiveAge)), 1.0), 0.0);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const {
        return sim::roundToTSFromDays( max(f1->sampleAgeOfDecay(rng), f2->sampleAgeOfDecay(rng)) );
    }

    unique_ptr<DecayFunction> hetSample(double hetFactor) const {
        unique_ptr<DecayFunction> f1hetSample = f1->hetSample(hetFactor);
        unique_ptr<DecayFunction> f2hetSample = f2->hetSample(hetFactor);
        unique_ptr<OperatorDecayFunction> copy = make_unique<OperatorDecayFunction>(*this, move(f1hetSample), move(f2hetSample));
        return move(copy);
    }

private:
    unique_ptr<DecayFunction> f1, f2;
    T op;
};

// -----  interface / static functions  -----

unique_ptr<DecayFunction> DecayFunction::makeObject(
    const scnXml::DecayFunction& elt, const char* eltName
){
    // Type mostly equivalent to a std::string:
    const scnXml::Function& func = elt.getFunction();
    if( func == "constant" ){
        return unique_ptr<DecayFunction>(new ConstantDecayFunction( elt ));
    }else if( func == "step" ){
        return unique_ptr<DecayFunction>(new StepDecayFunction( elt ));
    }else if( func == "linear" ){
        return unique_ptr<DecayFunction>(new LinearDecayFunction( elt ));
    }else if( func == "exponential" ){
        return unique_ptr<DecayFunction>(new ExponentialDecayFunction( elt ));
    }else if( func == "weibull" ){
        return unique_ptr<DecayFunction>(new WeibullDecayFunction( elt ));
    }else if( func == "hill" ){
        return unique_ptr<DecayFunction>(new HillDecayFunction( elt ));
    }else if( func == "smooth-compact" ){
        return unique_ptr<DecayFunction>(new SmoothCompactDecayFunction( elt ));
    }else if( func == "plus" ){
        return unique_ptr<DecayFunction>(new OperatorDecayFunction<std::plus<double>>( elt ));
    }else if( func == "minus" ){
        return unique_ptr<DecayFunction>(new OperatorDecayFunction<std::minus<double>>( elt ));
    }else if( func == "divides" ){
        return unique_ptr<DecayFunction>(new OperatorDecayFunction<std::divides<double>>( elt ));
    }else if( func == "multiplies" ){
        return unique_ptr<DecayFunction>(new OperatorDecayFunction<std::multiplies<double>>( elt ));
    }else{
        throw util::xml_scenario_error("decay function type " + string(func) + " of " + string(eltName) + " unrecognized");
    }
}

} }
