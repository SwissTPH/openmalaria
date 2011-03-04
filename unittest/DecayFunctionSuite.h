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

#ifndef Hmod_DecayFunctionSuite
#define Hmod_DecayFunctionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "util/DecayFunction.h"

using ::OM::util::DecayFunctionValue;
using ::OM::util::DecayFuncHet;

class DecayFunctionSuite : public CxxTest::TestSuite
{
public:
    DecayFunctionSuite () :
        dfElt( "", 10.0, 1.9 )
    {
        dfElt.setK( 1.6 );
        dHetOne.tMult = 1.0;
        dHetHalf.tMult = 0.5;
    }
    
    void setUp() {
        TimeStep::init(5,90.0);
    }
    
    void testBad() {
        dfElt.setFunction( "unknown" );
        string errMsg = "decay function type unknown of DecayFunctionSuite unrecognized";
        TS_ASSERT_THROWS_EQUALS( df.set( dfElt, "DecayFunctionSuite" ),
                          const std::runtime_error &e, e.what(), errMsg );
    }
    
    void testConstant () {
        dfElt.setFunction( "constant" );
        df.set( dfElt, "DecayFunctionSuite" );
        TS_ASSERT_APPROX( df.eval( TimeStep(0), dHetOne ), 1.9 );
        TS_ASSERT_APPROX( df.eval( TimeStep(1460), dHetHalf ), 1.9 );
    }
    
    void testStep () {
        dfElt.setFunction( "step" );
        df.set( dfElt, "DecayFunctionSuite" );
        TS_ASSERT_APPROX( df.eval( TimeStep(0), dHetOne ), 1.9 );
        TS_ASSERT_APPROX( df.eval( TimeStep(438), dHetOne ), 1.9 );
        TS_ASSERT_APPROX( df.eval( TimeStep(876), dHetHalf ), 1.9 );
        TS_ASSERT_APPROX( df.eval( TimeStep(1460), dHetOne ), 0.0 );
    }
    
    void testLinear () {
        dfElt.setFunction( "linear" );
        df.set( dfElt, "DecayFunctionSuite" );
        TS_ASSERT_APPROX( df.eval( TimeStep(0), dHetOne ), 1.9 );
        TS_ASSERT_APPROX( df.eval( TimeStep(438), dHetOne ), 0.76 );
        TS_ASSERT_APPROX( df.eval( TimeStep(876), dHetHalf ), 0.76 );
        TS_ASSERT_APPROX( df.eval( TimeStep(1460), dHetOne ), 0.0 );
    }
    
    void testExponential () {
        dfElt.setFunction( "exponential" );
        df.set( dfElt, "DecayFunctionSuite" );
        TS_ASSERT_APPROX( df.eval( TimeStep(0), dHetOne ), 1.9 );
        TS_ASSERT_APPROX( df.eval( TimeStep(438), dHetOne ), 1.2535325 );
        TS_ASSERT_APPROX( df.eval( TimeStep(876), dHetHalf ), 1.2535325 );
        TS_ASSERT_APPROX( df.eval( TimeStep(1460), dHetOne ), 0.475 );
    }
    
    void testWeibull () {
        dfElt.setFunction( "weibull" );
        df.set( dfElt, "DecayFunctionSuite" );
        TS_ASSERT_APPROX( df.eval( TimeStep(0), dHetOne ), 1.9 );
        TS_ASSERT_APPROX( df.eval( TimeStep(438), dHetOne ), 1.3989906 );
        TS_ASSERT_APPROX( df.eval( TimeStep(876), dHetHalf ), 1.3989906 );
        TS_ASSERT_APPROX( df.eval( TimeStep(1460), dHetOne ), 0.2323814 );
    }
    
    void testHill () {
        dfElt.setFunction( "hill" );
        df.set( dfElt, "DecayFunctionSuite" );
        TS_ASSERT_APPROX( df.eval( TimeStep(0), dHetOne ), 1.9 );
        TS_ASSERT_APPROX( df.eval( TimeStep(438), dHetOne ), 1.3179680 );
        TS_ASSERT_APPROX( df.eval( TimeStep(876), dHetHalf ), 1.3179680 );
        TS_ASSERT_APPROX( df.eval( TimeStep(1460), dHetOne ), 0.47129642 );
    }
    
    void testSmoothCompact () {
        dfElt.setFunction( "smooth-compact" );
        df.set( dfElt, "DecayFunctionSuite" );
        TS_ASSERT_APPROX( df.eval( TimeStep(0), dHetOne ), 1.9 );
        TS_ASSERT_APPROX( df.eval( TimeStep(438), dHetOne ), 0.77248235 );
        TS_ASSERT_APPROX( df.eval( TimeStep(876), dHetHalf ), 0.77248235 );
        TS_ASSERT_APPROX( df.eval( TimeStep(1460), dHetOne ), 0.0 );
    }
    
private:
    scnXml::DecayFunctionValue dfElt;
    DecayFunctionValue df;
    DecayFuncHet dHetOne;
    DecayFuncHet dHetHalf;
};

#endif
