/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "util/DecayFunction.h"
#include "util/errors.h"
#include "inputData.h"
#include "util/StreamValidator.h"
#include "util/random.h"

#include <cmath>
#include <stdexcept>
#include <boost/format.hpp>

namespace OM {
namespace util {

class BaseHetDecayFunction : public DecayFunction {
    double mu, sigma;
public:
    BaseHetDecayFunction( const scnXml::DecayFunction& elt ) :
        mu( elt.getMu() ), sigma( elt.getSigma() )
    {
        if( mu != 0.0 && sigma == 0.0 ){
            cerr << "Warning: for some decay function mu != 0 while sigma = 0 (or unspecified); in this case the mu parameter has no effect" << endl;
        }
    }
    
    virtual double getBaseTMult() const =0;
    
    DecayFuncHet hetSample () const{
        DecayFuncHet ret;
        if(sigma>0.0){
            ret.tMult = random::log_normal( mu, sigma );
        }else{
            assert(sigma==0.0);
            ret.tMult = 1.0;    // same answer as above but without using the random-number generator
        }
        ret.tMult *= getBaseTMult();
        return ret;
    }
    DecayFuncHet hetSample (NormalSample sample) const{
        DecayFuncHet ret;
        ret.tMult = sample.asLognormal( mu, sigma );
        ret.tMult *= getBaseTMult();
        return ret;
    }
};

class ConstantDecayFunction : public DecayFunction {
public:
    DecayFuncHet hetSample () const{
        return DecayFuncHet();
    }
    DecayFuncHet hetSample (NormalSample) const{
        return DecayFuncHet();
    }
    double eval(TimeStep age, DecayFuncHet sample) const{
        return 1.0;
    }
    TimeStep sampleAgeOfDecay () const{
        return TimeStep::never;
    }
};

class StepDecayFunction : public BaseHetDecayFunction {
public:
    StepDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        invL( 1.0 / (elt.getL() * TimeStep::stepsPerYear) )
    {}
    
    double getBaseTMult() const{
        return invL;
    }
    double eval(TimeStep age, DecayFuncHet sample) const{
        double effectiveAge = age.asInt() * sample.getTMult();
        if( effectiveAge < 1.0 ){
            return 1.0;
        }else{
            return 0.0;
        }
    }
    
    TimeStep sampleAgeOfDecay () const{
        return TimeStep::never;
    }
    
private:
    double invL;
};

class LinearDecayFunction : public BaseHetDecayFunction {
public:
    LinearDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        L( TimeStep::fromYears( elt.getL() ) ),
        invL( 1.0 / (elt.getL() * TimeStep::stepsPerYear) )
    {}
    
    double getBaseTMult() const{
        return invL;
    }
    double eval(TimeStep age, DecayFuncHet sample) const{
        double effectiveAge = age.asInt() * sample.getTMult();
        if( effectiveAge < 1.0 ){
            return 1.0 - effectiveAge;
        }else{
            return 0.0;
        }
    }
    
    TimeStep sampleAgeOfDecay () const{
        // Note: rounds to nearest. Object may decay instantly or at time L.
        return L * random::uniform_01();
    }
    
private:
    TimeStep L;
    double invL;
};

class ExponentialDecayFunction : public BaseHetDecayFunction {
public:
    ExponentialDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        negInvLambda( -log(2.0) / (elt.getL() * TimeStep::stepsPerYear) )
    {
        util::streamValidate(negInvLambda);
    }
    
    double getBaseTMult() const{
        return negInvLambda;
    }
    double eval(TimeStep age, DecayFuncHet sample) const{
        double effectiveAge = age.asInt() * sample.getTMult();
        return exp( effectiveAge );
    }
    
    TimeStep sampleAgeOfDecay () const{
        return TimeStep(log(random::uniform_01())/negInvLambda);
    }
    
private:
    double negInvLambda;
};

class WeibullDecayFunction : public BaseHetDecayFunction {
public:
    WeibullDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        constOverLambda( pow(log(2.0),1.0/elt.getK()) / (elt.getL() * TimeStep::stepsPerYear) ),
        k( elt.getK() )
    {}
    
    double getBaseTMult() const{
        return constOverLambda;
    }
    double eval(TimeStep age, DecayFuncHet sample) const{
        double effectiveAge = age.asInt() * sample.getTMult();
        return exp( -pow(effectiveAge, k) );
    }
    
    TimeStep sampleAgeOfDecay () const{
        return TimeStep( pow( -log(random::uniform_01()), 1.0/k ) / constOverLambda );
    }
    
private:
    double constOverLambda;
    double k;
};

class HillDecayFunction : public BaseHetDecayFunction {
public:
    HillDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        invL( 1.0 / (elt.getL() * TimeStep::stepsPerYear) ),
        k( elt.getK() )
    {}
    
    double getBaseTMult() const{
        return invL;
    }
    double eval(TimeStep age, DecayFuncHet sample) const{
        double effectiveAge = age.asInt() * sample.getTMult();
        return 1.0 / (1.0 + pow(effectiveAge, k));
    }
    
    TimeStep sampleAgeOfDecay () const{
        return TimeStep( pow( 1.0 / random::uniform_01() - 1.0, 1.0/k ) / invL );
    }
    
private:
    double invL, k;
};

class SmoothCompactDecayFunction : public BaseHetDecayFunction {
public:
    SmoothCompactDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        invL( 1.0 / (elt.getL() * TimeStep::stepsPerYear) ),
        k( elt.getK() )
    {}
    
    double getBaseTMult() const{
        return invL;
    }
    double eval(TimeStep age, DecayFuncHet sample) const{
        double effectiveAge = age.asInt() * sample.getTMult();
        if( effectiveAge < 1.0 ){
            return exp( k - k / (1.0 - pow(effectiveAge, 2.0)) );
        }else{
            return 0.0;
        }
    }
    
    TimeStep sampleAgeOfDecay () const{
        return TimeStep( sqrt( 1.0 - k / (k - log( random::uniform_01() )) ) / invL );
    }
    
private:
    double invL, k;
};


// -----  interface / static functions  -----

auto_ptr<DecayFunction> DecayFunction::makeObject(
    const scnXml::DecayFunction& elt, const char* eltName
){
    // Type mostly equivalent to a std::string:
    const scnXml::Function& func = elt.getFunction();
    if( func == "constant" ){
        return auto_ptr<DecayFunction>(new ConstantDecayFunction);
    }else if( func == "step" ){
        return auto_ptr<DecayFunction>(new StepDecayFunction( elt ));
    }else if( func == "linear" ){
        return auto_ptr<DecayFunction>(new LinearDecayFunction( elt ));
    }else if( func == "exponential" ){
        return auto_ptr<DecayFunction>(new ExponentialDecayFunction( elt ));
    }else if( func == "weibull" ){
        return auto_ptr<DecayFunction>(new WeibullDecayFunction( elt ));
    }else if( func == "hill" ){
        return auto_ptr<DecayFunction>(new HillDecayFunction( elt ));
    }else if( func == "smooth-compact" ){
        return auto_ptr<DecayFunction>(new SmoothCompactDecayFunction( elt ));
    }else{
        throw util::xml_scenario_error( (boost::format( "decay function type %1% of %2% unrecognized" ) %func %eltName).str() );
    }
}
auto_ptr<DecayFunction> DecayFunction::makeConstantObject(){
    return auto_ptr<DecayFunction>(new ConstantDecayFunction);
}


} }