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

#ifndef Hmod_DummyInfectionSuite
#define Hmod_DummyInfectionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "WithinHost/Infection/DummyInfection.h"
#include "WithinHost/CommonWithinHost.h"
#include <limits>

using namespace OM::WithinHost;

class DummyInfectionSuite : public CxxTest::TestSuite
{
public:
    void setUp () {
        UnittestUtil::initTime(1);
        UnittestUtil::Infection_init_latentP_and_NaN ();
        DummyInfection::init();
        // pkpdID (1st value) isn't important since we're not using drug model here:
        infection = CommonWithinHost::createInfection( 0xFFFFFFFF );
        for( SimTime d = sim::now(), end = sim::now() + sim::fromDays(15); d < end; d += sim::oneDay() ){
            // blood stage starts 15 days after creation
            UnittestUtil::incrTime( sim::oneDay() );
            infection->update( 1.0, d );
        }
    }
    void tearDown () {
        delete infection;
    }

    void testNewInf () {
        TS_ASSERT_APPROX (infection->getDensity(), 16.00000000288086086);
    }

    void testUpdatedInf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        TS_ASSERT_APPROX (infection->getDensity(), 128.00000008620828820);
    }
    void testUpdated2Inf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        TS_ASSERT_APPROX (infection->getDensity(), 1024.00000082264208600);
    }

    void testUpdatedReducedInf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (0.1, sim::now());
        // This is, as expected, 1/10th of that in testUpdated2Inf
        TS_ASSERT_APPROX (infection->getDensity(), 102.40000008226420860);
    }
    void testUpdatedReducedInf2 () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (0.1, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        // This is nearly the same
        TS_ASSERT_APPROX (infection->getDensity(), 102.00000008286288040);
    }
    
private:
    CommonInfection* infection;
};

#endif
