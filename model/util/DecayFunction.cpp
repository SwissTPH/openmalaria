/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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
        DecayFuncHet ret;
        ret.tMult = 1.0;
        return ret;
    }
    DecayFuncHet hetSample (NormalSample sample) const{
        return hetSample();
    }
    
    double eval(double effectiveAge) const{
        // Note: we now require all decay functions to return 0 when time > 0
        // and the DecayFuncHet is default-constructed. So const *after deployment*.
        if( effectiveAge == numeric_limits<double>::infinity() )
            return 0.0;
        return 1.0;
    }
    SimTime sampleAgeOfDecay () const{
        return sim::future();        // decay occurs "in the future" (don't use sim::never() because that is interpreted as being in the past)
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
    
    double getBaseTMult() const{
        return invL;
    }
    double eval(double effectiveAge) const{
        if( effectiveAge < 1.0 ){
            return 1.0;
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay () const{
        return sim::roundToTSFromDays( 1.0 / invL );
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
    
    double getBaseTMult() const{
        return invL;
    }
    double eval(double effectiveAge) const{
        if( effectiveAge < 1.0 ){
            return 1.0 - effectiveAge;
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay () const{
        // Note: rounds to nearest. Object may decay instantly or at time L.
        return sim::roundToTSFromDays(random::uniform_01() / invL);
    }
    
private:
    double invL;
};

class ExponentialDecayFunction : public BaseHetDecayFunction {
public:
    ExponentialDecayFunction( const scnXml::DecayFunction& elt ) :
        BaseHetDecayFunction( elt ),
        invLambda( log(2.0) / readLToDays(elt) )
    {
        util::streamValidate(invLambda);
    }
    
    double getBaseTMult() const{
        return invLambda;
    }
    double eval(double effectiveAge) const{
        return exp( -effectiveAge );
    }
    
    SimTime sampleAgeOfDecay () const{
        return sim::roundToTSFromDays( -log(random::uniform_01()) / invLambda );
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
    
    double getBaseTMult() const{
        return constOverLambda;
    }
    double eval(double effectiveAge) const{
        return exp( -pow(effectiveAge, k) );
    }
    
    SimTime sampleAgeOfDecay () const{
        return sim::roundToTSFromDays( pow( -log(random::uniform_01()), 1.0/k ) / constOverLambda );
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
    
    double getBaseTMult() const{
        return invL;
    }
    double eval(double effectiveAge) const{
        return 1.0 / (1.0 + pow(effectiveAge, k));
    }
    
    SimTime sampleAgeOfDecay () const{
        return sim::roundToTSFromDays( pow( 1.0 / random::uniform_01() - 1.0, 1.0/k ) / invL );
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
    
    double getBaseTMult() const{
        return invL;
    }
    double eval(double effectiveAge) const{
        if( effectiveAge < 1.0 ){
            return exp( k - k / (1.0 - pow(effectiveAge, 2.0)) );
        }else{
            return 0.0;
        }
    }
    
    SimTime sampleAgeOfDecay () const{
        return sim::roundToTSFromDays( sqrt( 1.0 - k / (k - log( random::uniform_01() )) ) / invL );
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