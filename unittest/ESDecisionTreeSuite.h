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

#ifndef Hmod_ESDecisionTreeSuite
#define Hmod_ESDecisionTreeSuite

#include <cxxtest/TestSuite.h>
#include "Clinical/ESCaseManagement.h"
#include "WithinHost/WithinHostModel.h"
#include "UnittestUtil.hpp"
#include <limits>

using namespace OM::Clinical;
using namespace OM::Pathogenesis;
using namespace OM::WithinHost;

class ESDecisionTreeSuite : public CxxTest::TestSuite
{
public:
    void setUp () {
	UnittestUtil::EmpiricalWHM_setup();
	// could seed random-number-generator, but shouldn't affect outcomes
	UnittestUtil::PkPdSuiteSetup (PkPd::PkPdModel::LSTM_PKPD);
	whm = WithinHostModel::createWithinHostModel();
    }
    void tearDown () {
	delete whm;
	UnittestUtil::PkPdSuiteTearDown ();
    }
    
    //TODO: test ESDecisionValue operations
    //TODO: test ESDecisionRandom
    
    void testUC2Test () {
	ESDecisionUC2Test d( dvMap );
	ESHostData hd(
		numeric_limits< double >::quiet_NaN(),
		*whm,
		STATE_MALARIA
	);
	TS_ASSERT_EQUALS( d.determine( ESDecisionValue(), hd ), dvMap.get( "case", "UC1" ) );
	hd.pgState = static_cast<State>( STATE_MALARIA | SECOND_CASE );
	TS_ASSERT_EQUALS( d.determine( ESDecisionValue(), hd ), dvMap.get( "case", "UC2" ) );
    }
    
    void testAge5Test () {
	ESDecisionAge5Test d( dvMap );
	ESHostData hd(
		4.99,
		*whm,
		STATE_MALARIA
	);
	TS_ASSERT_EQUALS( d.determine( ESDecisionValue(), hd ), dvMap.get( "age5Test", "under5" ) );
	hd.ageYears = 5.0;
	TS_ASSERT_EQUALS( d.determine( ESDecisionValue(), hd ), dvMap.get( "age5Test", "over5" ) );
    }
    
    void testParasiteTest () {
	ESDecisionParasiteTest d( dvMap );
	ESHostData hd(
		numeric_limits< double >::quiet_NaN(),
		*whm,
		STATE_MALARIA
	);
	//NOTE: the test isn't implemented yet, hence always returns negative
	TS_ASSERT_EQUALS( d.determine( dvMap.get( "test", "microscopy" ), hd ), dvMap.get( "result", "positive" ) );
	TS_ASSERT_EQUALS( d.determine( dvMap.get( "test", "RDT" ), hd ), dvMap.get( "result", "positive" ) );
	TS_ASSERT_EQUALS( d.determine( dvMap.get( "test", "none" ), hd ), ESDecisionValue() );
    }
    
private:
    ESDecisionValueMap dvMap;
    WithinHostModel* whm;
};

#endif
