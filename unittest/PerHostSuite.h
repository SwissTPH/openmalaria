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
// A very simple test for PerHostTransmission

#ifndef Hmod_PerHostSuite
#define Hmod_PerHostSuite

#include <cxxtest/TestSuite.h>
#include "Transmission/PerHost.h"

class PerHostSuite : public CxxTest::TestSuite
{
public:
  void setUp () {
    PerHostTransmission::initParameters();
  }
  
  void testRelativeAvailability () {
    TS_ASSERT_DELTA (PerHostTransmission::getRelativeAvailability(7.0), 0.51263046437755255, 0.00000000000000000);
  }
};

#endif
