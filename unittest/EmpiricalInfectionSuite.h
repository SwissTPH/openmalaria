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

#ifndef Hmod_EmpiricalInfectionSuite
#define Hmod_EmpiricalInfectionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "WithinHost/Infection/EmpiricalInfection.h"
#include "WithinHost/CommonWithinHost.h"
#include "util/random.h"
#include <limits>

using namespace OM::WithinHost;

class EmpiricalInfectionSuite : public CxxTest::TestSuite
{
public:
    void setUp () {
        UnittestUtil::initTime(1);
        UnittestUtil::Infection_init_latentP_and_NaN ();
        EmpiricalInfection::init();
        util::random::seed (83);	// seed is unimportant, but must be fixed
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
        TS_ASSERT_APPROX (infection->getDensity(), 0.00000000000000000);
    }

    // Parasite growth is stochastic, so there's not a lot we can test, except for reproducability
    void testUpdatedInf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        TS_ASSERT_APPROX (infection->getDensity(), 15.36758760023472284);
    }
    void testUpdated2Inf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        TS_ASSERT_APPROX (infection->getDensity(), 4.94261787639103382);
    }
    void testUpdated3Inf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        TS_ASSERT_APPROX (infection->getDensity(), 162.62062791268144860);
    }
    void testUpdated4Inf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        TS_ASSERT_APPROX (infection->getDensity(), 6.10393200785528424);
    }
    void testUpdatedInf1 () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        TS_ASSERT_APPROX (infection->getDensity(), 15.36758760023472284);
    }

    void testUpdatedReducedInf () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (0.1, sim::now());
        // This is, as expected, 1/10th of that in testUpdated2Inf
        TS_ASSERT_APPROX (infection->getDensity(), 0.49426178763910338);
    }
    void testUpdatedReducedInf2 () {
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (0.1, sim::now());
        UnittestUtil::incrTime( sim::oneTS() );
        infection->update (1.0, sim::now());
        // This is completely different due to stochasitic effects
        TS_ASSERT_APPROX (infection->getDensity(), 1.97582432565095644);
    }

private:
    CommonInfection* infection;
};

#endif
