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
// Unittest for the AgeGroupData

#ifndef Hmod_AgeGroupDataSuite
#define Hmod_AgeGroupDataSuite

#include <cxxtest/TestSuite.h>
#include "AgeGroupData.h"
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include <limits>

using namespace OM;

class AgeGroupDataSuite : public CxxTest::TestSuite
{
public:
    void setUp () {
	agd = new AgeGroupData;
    }
    void tearDown () {
	delete agd;
    }
    
    void testLowerbound() {
	double age = 0.0;
	agd->update( age );
	TS_ASSERT_APPROX (agd->ageToWeight(age), 13.9856718);
    }
    
    void testUpperbound() {
	double age = 100.0;
	agd->update( age );
	TS_ASSERT_APPROX (agd->ageToWeight(age), 60.0);
    }
    
    void testStandardCaseLowerbound(){
	double age = 14.99;
	agd->update( age );
	TS_ASSERT_APPROX (agd->ageToWeight(age), 49.48396092);
    }
    
    void testStandardCase() {
	double age = 17.0;
	agd->update( age );
	TS_ASSERT_APPROX (agd->ageToWeight(age), 51.44412863);
    }
    
    void testStandardCaseUpperbound(){
	double age = 19.99;
	agd->update( age );
	TS_ASSERT_APPROX (agd->ageToWeight(age), 54.36);
    }
    
    AgeGroupData *agd;
};

#endif
