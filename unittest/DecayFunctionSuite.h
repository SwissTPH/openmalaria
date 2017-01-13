/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
 
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

#ifndef Hmod_DecayFunctionSuite
#define Hmod_DecayFunctionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "util/DecayFunction.h"

using ::OM::util::DecayFunction;
using ::OM::util::DecayFuncHet;

class DecayFunctionSuite : public CxxTest::TestSuite
{
public:
    DecayFunctionSuite () :
        dfElt( "" )
    {
        dfElt.setL( "10y" );
        dfElt.setK( 1.6 );
    }
    
    void setUp() {
        UnittestUtil::initTime(5);
    }
    
    void testBad() {
        dfElt.setFunction( "unknown" );
        string errMsg = "decay function type unknown of DecayFunctionSuite unrecognized";
        TS_ASSERT_THROWS_EQUALS( df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" ),
                          const std::runtime_error &e, e.what(), errMsg );
    }
    
    void testConstant () {
        dfElt.setFunction( "constant" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        // First test: with a default-constructed DecayFuncHet and positive age, result should be zero
        DecayFuncHet dHet;
        TS_ASSERT_EQUALS( df->eval( SimTime::fromDays(5), dHet ), 0.0 );
        // Second test: with an appropriately sampled helper value, we should get the results we want
        dHet = df->hetSample();
        TS_ASSERT_APPROX( df->eval( SimTime::zero(), dHet ), 1.0 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(10), dHet ), 1.0 );
        // Third test: time of decay (age plus now) should always be in the future
        TS_ASSERT_EQUALS( df->sampleAgeOfDecay(), SimTime::future() );
    }
    
    void testStep () {
        dfElt.setFunction( "step" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        DecayFuncHet dHet;
        TS_ASSERT_EQUALS( df->eval( SimTime::fromDays(5), dHet ), 0.0 );
        dHet = df->hetSample();
        TS_ASSERT_APPROX( df->eval( SimTime::zero(), dHet ), 1.0 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(6), dHet ), 1.0 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(20), dHet ), 0.0 );
        TS_ASSERT_EQUALS( df->sampleAgeOfDecay(), SimTime::fromYearsI(10) );
    }
    
    void testLinear () {
        dfElt.setFunction( "linear" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        DecayFuncHet dHet;
        TS_ASSERT_EQUALS( df->eval( SimTime::fromDays(5), dHet ), 0.0 );
        dHet = df->hetSample();
        TS_ASSERT_APPROX( df->eval( SimTime::zero(), dHet ), 1.0 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(6), dHet ), 0.4 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(20), dHet ), 0.0 );
    }
    
    void testExponential () {
        dfElt.setFunction( "exponential" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        DecayFuncHet dHet;
        TS_ASSERT_EQUALS( df->eval( SimTime::fromDays(5), dHet ), 0.0 );
        dHet = df->hetSample();
        TS_ASSERT_APPROX( df->eval( SimTime::zero(), dHet ), 1.0 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(6), dHet ), 0.65975394736842108 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(20), dHet ), 0.25 );
    }
    
    void testWeibull () {
        dfElt.setFunction( "weibull" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        DecayFuncHet dHet = df->hetSample();
        TS_ASSERT_APPROX( df->eval( SimTime::zero(), dHet ), 1.0 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(6), dHet ), 0.73631084210526321 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(20), dHet ), 0.122306 );
    }
    
    void testHill () {
        dfElt.setFunction( "hill" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        DecayFuncHet dHet;
        TS_ASSERT_EQUALS( df->eval( SimTime::fromDays(5), dHet ), 0.0 );
        dHet = df->hetSample();
        TS_ASSERT_APPROX( df->eval( SimTime::zero(), dHet ), 1.0 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(6), dHet ), 0.6936673684210527 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(20), dHet ), 0.24805074736842106 );
    }
    
    void testSmoothCompact () {
        dfElt.setFunction( "smooth-compact" );
        df = DecayFunction::makeObject( dfElt, "DecayFunctionSuite" );
        DecayFuncHet dHet;
        TS_ASSERT_EQUALS( df->eval( SimTime::fromDays(5), dHet ), 0.0 );
        dHet = df->hetSample();
        TS_ASSERT_APPROX( df->eval( SimTime::zero(), dHet ), 1.0 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(6), dHet ), 0.40656965789473687 );
        TS_ASSERT_APPROX( df->eval( SimTime::fromYearsI(20), dHet ), 0.0 );
    }
    
private:
    scnXml::DecayFunction dfElt;
    auto_ptr<DecayFunction> df;
};

#endif
