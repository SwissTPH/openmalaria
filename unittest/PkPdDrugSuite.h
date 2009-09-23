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
// Unittest for the drug model

#ifndef Hmod_PkPdDrugSuite
#define Hmod_PkPdDrugSuite

#include <cxxtest/TestSuite.h>
#include "Drug/PkPdDrug.h"
#include "ExtraAsserts.h"

class PkPdDrugSuite : public CxxTest::TestSuite
{
public:
  PkPdDrugSuite () {
    Global::interval = 1;	// I think the drug model is always going to be used with an interval of 1 day.
    Global::modelVersion = INCLUDES_PK_PD;
    DrugModel::init ();
  }
  
  void setUp () {
    proxy = new PkPdDrug ();
    proteome = &ProteomeInstance::instances[0];
  }
  void tearDown () {
    delete proxy;
  }
  
  void testNone () {
    TS_ASSERT_EQUALS (proxy->getDrugFactor (proteome), 1.0);
  }
  
  void testCq () {
    proxy->medicate ("CQ", 250000, 0, 21, 60);
    TS_ASSERT_APPROX (proxy->getDrugFactor (proteome), 0.12794177390385896);
  }
  
  void testCqDecayed () {
    proxy->medicate ("CQ", 250000, 0, 21, 60);
    proxy->decayDrugs ();
    TS_ASSERT_APPROX (proxy->getDrugFactor (proteome), 0.13760869542580346);
  }
  
  void testCq2Doses () {
    proxy->medicate ("CQ", 250000, 0, 60);
    proxy->decayDrugs ();
    proxy->medicate ("CQ", 250000, 0, 60);
    TS_ASSERT_APPROX (proxy->getDrugFactor (proteome), 0.07150144786339767);
  }
  
  PkPdDrug *proxy;
  ProteomeInstance *proteome;
};

#endif
