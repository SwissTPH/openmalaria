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
#include <stdexcept>

namespace OM {
namespace util {

/** Default value: should make all eval() calls return 0 (i.e. infinitely
 * old deployment). */
DecayFunctionHet::DecayFunctionHet() :  timeFactorHet(numeric_limits<double>::infinity()), sample(new NullDecayFunction()) {}



class BaseHetDecayFunction : public DecayFunction {
    LognormalSampler het;
public:
    BaseHetDecayFunction( const scnXml::DecayFunction& elt ) :
        DecayFunction( elt )
    {
        het.setMeanCV( 1.0, elt.getCV() );
    }
    
    virtual double getBaseTimeFactorHet() const =0;
    
    virtual DecayFunctionHet _sample(double s) const = 0;

    DecayFunctionHet hetSample (LocalRng& rng) const{
        //return DecayFunctionHet(het.sample(rng) * getBaseTimeFactorHet());
        return _sample(het.sample(rng) * getBaseTimeFactorHet());
    }
    DecayFunctionHet hetSample (NormalSample sample) const{
        //return DecayFunctionHet(het.sample(sample) * getBaseTimeFactorHet());
        return _sample(het.sample(sample) * getBaseTimeFactorHet());
    }
};

class ConstantDecayFunction : public BaseHetDecayFunction {
public:
    ConstantDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt )
    {}

    double getBaseTimeFactorHet() const{
        return 1.0;
    }
    
    double eval(double effectiveAge) const{
        // Note: we now require all decay functions to return 0 when time > 0
        // and the DecayFunctionHet is default-constructed. So const *after deployment*.
        if( effectiveAge == numeric_limits<double>::infinity() )
            return 0.0;
        return 1.0;
    }
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::future();        // decay occurs "in the future" (don't use sim::never() because that is interpreted as being in the past)
    }

    DecayFunctionHet _sample(double s) const {
        shared_ptr<DecayFunction> copy = make_shared<ConstantDecayFunction>(*this);
        return DecayFunctionHet(s, copy);
    }
};

double readLToDays( const scnXml::DecayFunction& elt ){
    if( !elt.getL().present() ){
        throw util::xml_scenario_error( "decay function: attribute L required" );
    }
    return UnitParse::durationToDays(elt.getL().get(), UnitParse::YEARS);
}

class StepDecayFunction : public BaseHetDecayFunction {
public:
    StepDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        invL( 1.0 / readLToDays(elt) )
    {}
    
    double getBaseTimeFactorHet() const{
        return invL;
    }
    double eval(double effectiveAge) const{
        if( effectiveAge < 1.0 ){
            return 1.0;
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( 1.0 / invL );
    }

    DecayFunctionHet _sample(double s) const {
        shared_ptr<DecayFunction> copy = make_shared<StepDecayFunction>(*this);
                cout << "creating LinearDecayFunction" << endl;
        return DecayFunctionHet(s, copy);
    }
    
private:
    double invL;        // 1 / days
};

class LinearDecayFunction : public BaseHetDecayFunction {
public:
    LinearDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        invL( 1.0 / readLToDays(elt) )
    {}
    
    double getBaseTimeFactorHet() const{
        return invL;
    }
    double eval(double effectiveAge) const{
        if( effectiveAge < 1.0 ){
            return 1.0 - effectiveAge;
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        // Note: rounds to nearest. Object may decay instantly or at time L.
        return sim::roundToTSFromDays(rng.uniform_01() / invL);
    }

    DecayFunctionHet _sample(double s) const {
        shared_ptr<DecayFunction> copy = make_shared<LinearDecayFunction>(*this);
        return DecayFunctionHet(s, copy);
    }
    
private:
    double invL;
};

class ExponentialDecayFunction : public BaseHetDecayFunction {
public:
    ExponentialDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        invLambda( log(2.0) / readLToDays(elt) )
    {}
    
    double getBaseTimeFactorHet() const{
        return invLambda;
    }
    double eval(double effectiveAge) const{
        return exp( -effectiveAge );
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( -log(rng.uniform_01()) / invLambda );
    }

    DecayFunctionHet _sample(double s) const {
        shared_ptr<DecayFunction> copy = make_shared<ExponentialDecayFunction>(*this);
        return DecayFunctionHet(s, copy);
    }
    
private:
    double invLambda;
};

class WeibullDecayFunction : public BaseHetDecayFunction {
public:
    WeibullDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        constOverLambda( pow(log(2.0),1.0/elt.getK()) / readLToDays(elt) ),
        k( elt.getK() )
    {}
    
    double getBaseTimeFactorHet() const{
        return constOverLambda;
    }
    double eval(double effectiveAge) const{
        double p = -pow(effectiveAge, k);
        if(p < -700.0)
            return 0.0;
        else
            return exp(p);
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( pow( -log(rng.uniform_01()), 1.0/k ) / constOverLambda );
    }

    DecayFunctionHet _sample(double s) const {
        shared_ptr<DecayFunction> copy = make_shared<WeibullDecayFunction>(*this);
        return DecayFunctionHet(s, copy);
    }
    
private:
    double constOverLambda;
    double k;
};

class HillDecayFunction : public BaseHetDecayFunction {
public:
    HillDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        invL( 1.0 / readLToDays(elt) ),
        k( elt.getK() )
    {}
    
    double getBaseTimeFactorHet() const{
        return invL;
    }
    double eval(double effectiveAge) const{
        return 1.0 / (1.0 + pow(effectiveAge, k));
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( pow( 1.0 / rng.uniform_01() - 1.0, 1.0/k ) / invL );
    }

    DecayFunctionHet _sample(double s) const {
        shared_ptr<DecayFunction> copy = make_shared<HillDecayFunction>(*this);
        return DecayFunctionHet(s, copy);
    }
    
private:
    double invL, k;
};

class SmoothCompactDecayFunction : public BaseHetDecayFunction {
public:
    SmoothCompactDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        invL( 1.0 / readLToDays(elt) ),
        k( elt.getK() )
    {}
    
    double getBaseTimeFactorHet() const{
        return invL;
    }
    double eval(double effectiveAge) const{
        if( effectiveAge < 1.0 ){
            return exp( k - k / (1.0 - pow(effectiveAge, 2.0)) );
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay (LocalRng& rng) const{
        return sim::roundToTSFromDays( sqrt( 1.0 - k / (k - log( rng.uniform_01() )) ) / invL );
    }

    DecayFunctionHet _sample(double s) const {
        shared_ptr<DecayFunction> copy = make_shared<SmoothCompactDecayFunction>(*this);
        return DecayFunctionHet(s, copy);
    }
    
private:
    double invL, k;
};

// class BiphasicDecayFunction : public BaseHetDecayFunction<BiphasicDecayFunction> {
// public:
//     BiphasicDecayFunction( const scnXml::DecayFunction& elt ) :
//         BaseHetDecayFunction( elt )
//     {
//         if( !elt.getInitial_efficacy().present() )
//             throw util::xml_scenario_error( "biphasic decay function: attribute initial_efficacy required" );
//         if( !elt.getRho().present() )
//             throw util::xml_scenario_error( "biphasic decay function: attribute rho required" );
//         if( !elt.getHalflife_long().present() )
//             throw util::xml_scenario_error( "biphasic decay function: attribute halflife_long required" );
//         if( !elt.getHalflife_short().present() )
//             throw util::xml_scenario_error( "biphasic decay function: attribute halflife_short required" );

//         initial_efficacy = elt.getInitial_efficacy().get();
//         rho = elt.getRho().get();
//         halflife_long = UnitParse::durationToDays(elt.getHalflife_long().get(), UnitParse::YEARS);
//         halflife_short = UnitParse::durationToDays(elt.getHalflife_short().get(), UnitParse::YEARS);

//         invL1 = log(2.0) / halflife_short;
//         invL2 = log(2.0) / halflife_long;
//     }
    
//     double getBaseTimeFactorHet() const{
//         return 1.0;
//     }
//     double eval(double effectiveAge) const{
//         return initial_efficacy * (
//             rho * exp(-effectiveAge * invL1) 
//             + (1.0 - rho) * exp(-effectiveAge * invL2)
//         );
//     }
    
//     SimTime sampleAgeOfDecay (LocalRng& rng) const{
//         throw util::xml_scenario_error( "biphasic decay function: sampleAgeOfDecay not implemented" );
//     }
    
// private:
//     double invL1, invL2;
//     double initial_efficacy;
//     double rho;
//     double halflife_long;
//     double halflife_short;
// };

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
    // }else if( func == "biphasic" ){
    //     return unique_ptr<DecayFunction>(new BiphasicDecayFunction( elt ));
    }else{
        throw util::xml_scenario_error("decay function type " + string(func) + " of " + string(eltName) + " unrecognized");
    }
}

} }
