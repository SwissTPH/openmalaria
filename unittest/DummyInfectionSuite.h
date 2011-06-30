/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#ifndef Hmod_DummyInfectionSuite
#define Hmod_DummyInfectionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "WithinHost/Infection/DummyInfection.h"
#include <limits>

using namespace OM::WithinHost;

class DummyInfectionSuite : public CxxTest::TestSuite
{
public:
    void setUp () {
	UnittestUtil::Infection_init_NaN ();
	DummyInfection::init();
	TimeStep::init( 1, 90.0 );
	TimeStep::simulation = TimeStep(1);	// value isn't really important
	infection = new DummyInfection (TimeStep::simulation, 0xFFFFFFFF);	// pkpdID (1st value) isn't important since we're not using drug model here
        TimeStep::simulation = TimeStep(2);     // next time step
    }
    void tearDown () {
	delete infection;
    }
    
    void testNewInf () {
	TS_ASSERT_APPROX (infection->getDensity(), 16.00000000288086086);
    }
    
    void testUpdatedInf () {
	infection->update (1.0);
	TS_ASSERT_APPROX (infection->getDensity(), 128.00000008620828820);
    }
    void testUpdated2Inf () {
	infection->update (1.0);
        TimeStep::simulation += TimeStep(1);
	infection->update (1.0);
	TS_ASSERT_APPROX (infection->getDensity(), 1024.00000082264208600);
    }
    
    void testUpdatedReducedInf () {
	infection->update (1.0);
        TimeStep::simulation += TimeStep(1);
	infection->update (0.1);
	// This is, as expected, 1/10th of that in testUpdated2Inf
	TS_ASSERT_APPROX (infection->getDensity(), 102.40000008226420860);
    }
    void testUpdatedReducedInf2 () {
	infection->update (0.1);
        TimeStep::simulation += TimeStep(1);
	infection->update (1.0);
	// This is nearly the same
	TS_ASSERT_APPROX (infection->getDensity(), 102.00000008286288040);
    }
    
private:
    CommonInfection* infection;
};

#endif
