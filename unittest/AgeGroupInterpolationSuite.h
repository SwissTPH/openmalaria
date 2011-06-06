/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#ifndef Hmod_AgeGroupInterpolationSuite
#define Hmod_AgeGroupInterpolationSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "util/AgeGroupInterpolation.h"

using util::AgeGroupInterpolation;

class AgeGroupInterpolationSuite : public CxxTest::TestSuite
{
public:
    AgeGroupInterpolationSuite () {
        // Set maximum age to 90 years:
        UnittestUtil::AgeGroupInterpolation_init();
        agvElt = auto_ptr<scnXml::AgeGroupValues>(new scnXml::AgeGroupValues());
        scnXml::AgeGroupValues::GroupSequence& seq = agvElt->getGroup();
        seq.resize( dataLen, scnXml::Group(0.0,0.0) );
        for( size_t i = 0; i < dataLen; ++i ){
            seq[ i ].setLowerbound( stdLbounds[ i ] );
            seq[ i ].setValue( stdValues[ i ] );
        }
    }
    ~AgeGroupInterpolationSuite () {
    }
    
    void setUp () {
    }
    void tearDown () {
    }
    
    void testDummy() {
        AgeGroupInterpolation *o = AgeGroupInterpolation::dummyObject();
        TS_ASSERT_THROWS( o->eval(5.7), const std::logic_error &e );
        AgeGroupInterpolation::freeObject(o);
    }
    
    void testPiecewiseConst () {
        agvElt->setInterpolation( "none" );
        AgeGroupInterpolation *o = AgeGroupInterpolation::makeObject( *agvElt, "testPiecewiseConst" );
        for( size_t i = 0; i < testLen; ++i ){
            TS_ASSERT_APPROX( o->eval( testAges[ i ] ), piecewiseConstValues[ i ] );
        }
        AgeGroupInterpolation::freeObject(o);
    }
    
    void testLinearInterp () {
        agvElt->setInterpolation( "linear" );
        AgeGroupInterpolation *o = AgeGroupInterpolation::makeObject( *agvElt, "testLinearInterp" );
        for( size_t i = 0; i < testLen; ++i ){
            TS_ASSERT_APPROX( o->eval( testAges[ i ] ), linearInterpValues[ i ] );
        }
        AgeGroupInterpolation::freeObject(o);
    }
    
private:
    static const size_t dataLen = 5;
    static const size_t testLen = 8;
    static const double stdLbounds[dataLen];
    static const double stdValues[dataLen];
    static const double testAges[testLen];
    static const double piecewiseConstValues[testLen];
    static const double linearInterpValues[testLen];
    auto_ptr<scnXml::AgeGroupValues> agvElt;
};

const double AgeGroupInterpolationSuite::stdLbounds[] = {
    0,5,10,15,60
};
const double AgeGroupInterpolationSuite::stdValues[] = {
    6.08,3.81,2.62,4.05,5.41
};
const double AgeGroupInterpolationSuite::testAges[] = {
    // various ages, designed to test limits, boundary points and interpolation
    15.2,18.09,7.0,2.5,
    0.0,20.0,900.0,62.0
};
const double AgeGroupInterpolationSuite::piecewiseConstValues[] = {
    4.0499999999999998, 4.0499999999999998, 3.8100000000000001, 6.0800000000000001,
    6.0800000000000001, 4.0499999999999998, 5.4100000000000001, 5.4100000000000001
};
const double AgeGroupInterpolationSuite::linearInterpValues[] = {
    2.7744400000000002, 2.9397479999999998, 4.0369999999999999, 6.0800000000000001,
    6.0800000000000001, 3.0489999999999999, 5.4100000000000001, 4.938533333333333
};

#endif
