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

#include <cmath>
#include <stdexcept>
#include <boost/format.hpp>

namespace OM {
namespace util {

class ConstantDecayFunction : public DecayFunction {
public:
    double eval(TimeStep age) const{
        return 1.0;
    }
};

class LinearDecayFunction : public DecayFunction {
public:
    LinearDecayFunction( const scnXml::DecayFunction& elt ) :
        L( TimeStep::fromYears( elt.getL() ) ),
        invL( 1.0 / (elt.getL() * TimeStep::stepsPerYear) )
    {}
    
    double eval(TimeStep age) const{
        if( age < L ){
            return 1.0 - age.asInt() * invL;
        }else{
            return 0.0;
        }
    }
    
private:
    TimeStep L;
    double invL;
};

class ExponentialDecayFunction : public DecayFunction {
public:
    ExponentialDecayFunction( const scnXml::DecayFunction& elt ) :
        negInvLambda( -log(2.0) / (elt.getL() * TimeStep::stepsPerYear) )
    {
        util::streamValidate(negInvLambda);
    }
    
    double eval(TimeStep age) const{
        return exp( age.asInt() * negInvLambda );
    }
    
private:
    double negInvLambda;
};

class WeibullDecayFunction : public DecayFunction {
public:
    WeibullDecayFunction( const scnXml::DecayFunction& elt ) :
        constOverLambda( pow(log(2.0),1.0/elt.getK()) / (elt.getL() * TimeStep::stepsPerYear) ),
        k( elt.getK() )
    {}
    
    double eval(TimeStep age) const{
        return exp( -pow(age.asInt() * constOverLambda, k) );
    }
    
private:
    double constOverLambda;
    double k;
};

class HillDecayFunction : public DecayFunction {
public:
    HillDecayFunction( const scnXml::DecayFunction& elt ) :
        invL( 1.0 / (elt.getL() * TimeStep::stepsPerYear) ),
        k( elt.getK() )
    {}
    
    double eval(TimeStep age) const{
        return 1.0 / (1.0 + pow(age.asInt() * invL, k));
    }
    
private:
    double invL, k;
};

class ChitnisDecayFunction : public DecayFunction {
public:
    ChitnisDecayFunction( const scnXml::DecayFunction& elt ) :
        L( TimeStep::fromYears( elt.getL() ) ),
        invL( 1.0 / (elt.getL() * TimeStep::stepsPerYear) ),
        k( elt.getK() )
    {}
    
    double eval(TimeStep age) const{
        if( age < L ){
            return exp( k - k / (1.0 - pow(age.asInt() * invL, 2.0)) );
        }else{
            return 0.0;
        }
    }
    
private:
    TimeStep L;
    double invL, k;
};


// -----  interface / static functions  -----

shared_ptr<DecayFunction> DecayFunction::makeObject(
    const scnXml::DecayFunction& elt, const char* eltName
){
    // Type mostly equivalent to a std::string:
    const scnXml::Function& func = elt.getFunction();
    if( func == "constant" ){
        return shared_ptr<DecayFunction>(new ConstantDecayFunction);
    }else if( func == "linear" ){
        return shared_ptr<DecayFunction>(new LinearDecayFunction( elt ));
    }else if( func == "exponential" ){
        return shared_ptr<DecayFunction>(new ExponentialDecayFunction( elt ));
    }else if( func == "weibull" ){
        return shared_ptr<DecayFunction>(new WeibullDecayFunction( elt ));
    }else if( func == "hill" ){
        return shared_ptr<DecayFunction>(new HillDecayFunction( elt ));
    }else if( func == "chitnis" ){
        return shared_ptr<DecayFunction>(new ChitnisDecayFunction( elt ));
    }else{
        throw util::xml_scenario_error( (boost::format( "decay function type %1% of %2% unrecognized" ) %func %eltName).str() );
    }
}
shared_ptr<DecayFunction> DecayFunction::makeConstantObject(){
    return shared_ptr<DecayFunction>(new ConstantDecayFunction);
}


// -----  DecayFunctionValue  -----

void DecayFunctionValue::set (const scnXml::DecayFunctionValue& elt, const char* eltName){
    initial = elt.getInitial();
    decayFunc = DecayFunction::makeObject( elt, eltName );
}


} }