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
// Unittest for the EventScheduler case management

#ifndef Hmod_ESCaseManagementSuite
#define Hmod_ESCaseManagementSuite

#include <cxxtest/TestSuite.h>
#include "Clinical/ESCaseManagement.h"
#include "ExtraAsserts.h"
#include <limits>

using namespace OM::Clinical;

class ESCaseManagementSuite : public CxxTest::TestSuite
{
public:
    //TODO: test tree execution (that all necessary decisions are evaluated and outputs conglomerated as expected)
    
    //TODO: test trees handle "void" output correctly
    
    //TODO: test all decisions of ESDecisionTree.h
    
    //TODO: test treatments get matched from inputs
    
    //TODO: test treatment modifiers get applied and matched correctly: delay, qty multiplier, adherence
};

#endif
