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

#ifndef Hmod_InfectionImmunitySuite
#define Hmod_InfectionImmunitySuite

#include <cxxtest/TestSuite.h>
#include "WithinHost/DummyInfection.h"
#include <limits>

using namespace OM::WithinHost;

class InfectionImmunitySuite : public CxxTest::TestSuite
{
public:
    InfectionImmunitySuite () {
	UnittestUtil::Infection_init ();
    }
    
    void setUp () {
	infection = new DummyInfection (0xFFFFFFFF);	// pkpdID (1st value) isn't important since we're not using drug model here
    }
    void tearDown () {
	delete infection;
    }
    
    void testImmunity () {
	// Base case: virtually no immunity due to mother's immunity or past infections
	TS_ASSERT_APPROX (infection->immunitySurvivalFactor (100., 0., 0.), 1.0);
	// maternal immunity
	TS_ASSERT_APPROX (infection->immunitySurvivalFactor (0.1, 0., 0.), 0.30631938551881299);
	// past infections, no density
	TS_ASSERT_APPROX (infection->immunitySurvivalFactor (100., 100., 0.), 0.41995608739475931);
	// previous with cum. density of 1e8 (but none from current infection)
	TS_ASSERT_APPROX (infection->immunitySurvivalFactor (100., 100., 1e8), 0.17081918453312689);
    }
    
private:
    DummyInfection* infection;
};

#endif
