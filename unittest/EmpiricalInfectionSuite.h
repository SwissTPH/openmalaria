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

#ifndef Hmod_EmpiricalInfectionSuite
#define Hmod_EmpiricalInfectionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.hpp"
#include "WithinHost/EmpiricalInfection.h"
#include "util/gsl.h"
#include <limits>

using namespace OM::WithinHost;

class EmpiricalInfectionSuite : public CxxTest::TestSuite
{
public:
    EmpiricalInfectionSuite () {
	UnittestUtil::Infection_init_NaN ();
	EmpiricalInfection::initParameters();
    }
    
    void setUp () {
	gsl::setUp (83);	// seed is unimportant, but must be fixed
	Global::simulationTime = 1;	// value isn't really important
	infection = new EmpiricalInfection (0xFFFFFFFF, 1);	// pkpdID (1st value) isn't important since we're not using drug model here
    }
    void tearDown () {
	delete infection;
	gsl::tearDown ();
    }
    
    void testNewInf () {
	TS_ASSERT_APPROX (infection->getDensity(), 0.00000000000000000);
    }
    
    // Parasite growth is stochastic, so there's not a lot we can test, except for reproducability
    void testUpdatedInf () {
	infection->updateDensity (Global::simulationTime + 1, 1.0);
	TS_ASSERT_APPROX (infection->getDensity(), 15.36758760023472284);
    }
    void testUpdated2Inf () {
	infection->updateDensity (Global::simulationTime + 1, 1.0);
	infection->updateDensity (Global::simulationTime + 2, 1.0);
	TS_ASSERT_APPROX (infection->getDensity(), 4.94261787639103382);
    }
    void testUpdated3Inf () {
	infection->updateDensity (Global::simulationTime + 1, 1.0);
	infection->updateDensity (Global::simulationTime + 2, 1.0);
	infection->updateDensity (Global::simulationTime + 3, 1.0);
	TS_ASSERT_APPROX (infection->getDensity(), 162.62062791268144860);
    }
    void testUpdated4Inf () {
	infection->updateDensity (Global::simulationTime + 1, 1.0);
	infection->updateDensity (Global::simulationTime + 2, 1.0);
	infection->updateDensity (Global::simulationTime + 3, 1.0);
	infection->updateDensity (Global::simulationTime + 4, 1.0);
	TS_ASSERT_APPROX (infection->getDensity(), 6.10393200785528424);
    }
    void testUpdatedInf1 () {
	infection->updateDensity (Global::simulationTime + 1, 1.0);
	TS_ASSERT_APPROX (infection->getDensity(), 15.36758760023472284);
    }
    
    void testUpdatedReducedInf () {
	infection->updateDensity (Global::simulationTime + 1, 1.0);
	infection->updateDensity (Global::simulationTime + 2, 0.1);
	// This is, as expected, 1/10th of that in testUpdated2Inf
	TS_ASSERT_APPROX (infection->getDensity(), 0.49426178763910338);
    }
    void testUpdatedReducedInf2 () {
	infection->updateDensity (Global::simulationTime + 1, 0.1);
	infection->updateDensity (Global::simulationTime + 2, 1.0);
	// This is completely different due to stochasitic effects
	TS_ASSERT_APPROX (infection->getDensity(), 1.97582432565095644);
    }
    
private:
    CommonInfection* infection;
};

#endif
