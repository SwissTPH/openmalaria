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

#ifndef Hmod_HoshenPkPdSuite
#define Hmod_HoshenPkPdSuite

#include <cxxtest/TestSuite.h>
#include "PkPd/HoshenPkPdModel.h"
#include "ExtraAsserts.h"
#include "Global.h"
#include "util/ModelOptions.hpp"
#include <limits>

using namespace OM;
using namespace OM::PkPd;

class HoshenPkPdSuite : public CxxTest::TestSuite
{
public:
  HoshenPkPdSuite () {
    Global::interval = 1;	// I think the drug model is always going to be used with an interval of 1 day.
    util::ModelOptions::optSet.set (util::INCLUDES_PK_PD);
    
    //Note: we fudge this call since it's not so easy to falsely initialize scenario element.
    ProteomeManager::init ();
    //PkPdModel::init ();
    PkPdModel::activeModel = PkPdModel::HOSHEN_PKPD;
    HoshenDrugType::init();
    HoshenDrug::init ();
    HoshenPkPdModel::init();
  }
  
  void setUp () {
    proxy = new HoshenPkPdModel ();
    proteome_ID = ProteomeInstance::getInstances()[0].getProteomeID();
  }
  void tearDown () {
    delete proxy;
  }
  
  void testNone () {
      // Note: second parameter is age but shouldn't be used; to check this pass an NaN
    TS_ASSERT_EQUALS (proxy->getDrugFactor (proteome_ID, std::numeric_limits< double >::quiet_NaN()), 1.0);
  }
  
  void testCq () {
    proxy->medicate ("CQ", 250000, 0, 21);
    TS_ASSERT_APPROX (proxy->getDrugFactor (proteome_ID, std::numeric_limits< double >::quiet_NaN()), 0.12427429993973554);
  }
  
  void testCqDecayed () {
    proxy->medicate ("CQ", 250000, 0, 21);
    proxy->decayDrugs (21);
    TS_ASSERT_APPROX (proxy->getDrugFactor (proteome_ID, std::numeric_limits< double >::quiet_NaN()), 0.12608995630400068);
  }
  
  void testCq2Doses () {
    proxy->medicate ("CQ", 250000, 0, 21);
    proxy->decayDrugs (21);
    proxy->medicate ("CQ", 250000, 0, 21);
    TS_ASSERT_APPROX (proxy->getDrugFactor (proteome_ID, std::numeric_limits< double >::quiet_NaN()), 0.06809903879225410);
  }
  
  HoshenPkPdModel *proxy;
  uint32_t proteome_ID;
};

#endif
