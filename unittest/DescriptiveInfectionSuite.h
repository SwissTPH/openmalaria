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

#ifndef Hmod_DescriptiveInfectionSuite
#define Hmod_DescriptiveInfectionSuite

#include <cxxtest/TestSuite.h>
#include "WithinHost/DescriptiveInfection.h"
#include "util/gsl.h"
#include <limits>

using namespace OM::WithinHost;

class DescriptiveInfectionSuite : public CxxTest::TestSuite
{
public:
    DescriptiveInfectionSuite () {
	UnittestUtil::Infection_init ();
	UnittestUtil::DescriptiveInfection_init ();
	DescriptiveInfection::initParameters();	//FIXME: this gets data from InputData
    }
    
    void setUp () {
	gsl::setUp (83);	// seed is unimportant, but must be fixed
	Global::interval = 5;
	Global::simulationTime = 1;	// value isn't really important
	infection = new DescriptiveInfection ();
    }
    void tearDown () {
	delete infection;
	gsl::tearDown ();
    }
    
    // Should test DescriptiveInfection::determineDensities (and maybe others), but due to large number of inputs & 2 outputs isn't so easy.
    
private:
    DescriptiveInfection* infection;
};

#endif
