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
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

using namespace OM::Clinical;
using namespace OM::Pathogenesis;
using namespace OM::WithinHost;
using namespace boost::assign; // bring 'operator+=()' into scope

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
    
    void testESDecisionValue () {
	ESDecisionValueMap dvMap;	// use a separate map for this test
	
	vector<string> vals;
	vals += "1","2";
	dvMap.add_decision_values( "a", vals );
	
	// Assert getting unknown values/decisions fails:
	TS_ASSERT_THROWS_EQUALS( dvMap.get( "b", "1" ), const std::runtime_error &e, string(e.what()), "ESDecisionValueMap::get(): no decision b" );
	TS_ASSERT_THROWS_EQUALS( dvMap.get( "a", "4" ), const std::runtime_error &e, string(e.what()), "ESDecisionValueMap::get(): no value a(4)" );
	
	// check we can re-add the same decision:
	TS_ASSERT_THROWS_NOTHING( dvMap.add_decision_values( "a", vals ) );
	// but not if it's different:
	vals[1] = "5";
	TS_ASSERT_THROWS_EQUALS( dvMap.add_decision_values( "a", vals ), const std::runtime_error &e, string(e.what()), "CaseManagement: decision a's values don't match; expected value: 2" );
	vals.resize(4,"2");
	vals[3]="6";	// now "1","5","2","6"
	TS_ASSERT_THROWS_EQUALS( dvMap.add_decision_values( "a", vals ), const std::runtime_error &e, string(e.what()), "CaseManagement: decision a's values don't match; unexpected values: 5 6" );
	
	// check we can't use void:
	vals.resize(1);
	vals[0] = "void";
	TS_ASSERT_THROWS_EQUALS( dvMap.add_decision_values( "v", vals ), const std::runtime_error &e, string(e.what()), "void can not be a declared output of a decision" );
	
	vals.resize(2,"9");
	vals[0] = "9";	// "9","9"
	TS_ASSERT_THROWS_EQUALS( dvMap.add_decision_values( "c", vals ), const std::runtime_error &e, string(e.what()), "CaseManagement: decision c's value 9 in value list twice!" );
	
	vals.resize(1);
	dvMap.add_decision_values( "c", vals );
	
	// check what outcomes we get:
	typedef ESDecisionValue::id_type id_type;
	TS_ASSERT_EQUALS( dvMap.get( "a", "1" ).id, static_cast<id_type>(1) );
	TS_ASSERT_EQUALS( dvMap.get( "a", "2" ).id, static_cast<id_type>(2) );
	TS_ASSERT_EQUALS( dvMap.get( "c", "9" ).id, static_cast<id_type>(4) );
	
	ostringstream val_string;
	val_string << dvMap.format( dvMap.get( "a", "1" ) | dvMap.get( "c", "9" ) );
	TS_ASSERT_EQUALS( val_string.str(), "a(1), c(9)" );
    }
    
    /* Runs d.determine( input, hd ) N times
     * returns the proportion of these runs where the output equalled expectedOutput
     */
    double determineNTimes (int N, const ESDecisionTree& d, const ESDecisionValue input, const ESHostData& hd, const ESDecisionValue expectedOutput) {
	int nExpected = 0;
	for (int i = 0; i < N; ++i) {
	    if( d.determine( input, hd ) == expectedOutput )
		++nExpected;
	}
	return double(nExpected) / double(N);
    }
    
    void testRandom () {
	ESHostData hd(
		numeric_limits< double >::quiet_NaN(),
		*whm,
		NONE
	);
	
	// random decision
	scnXml::HSESDecision ut_r_xml ("\
		p(.5) {\
		    p(.9): a\
		    p(.1): b\
		}\
		p(.5) {\
		    p(.7): a\
		    p(.3): b\
		}",
	    "myR",	// decision
	    "",	// depends
	    "a,b"	// values
	);
	ESDecisionRandom ut_r( dvMap, ut_r_xml );
	
	const int N = 10000;
	const double LIM = .02;
	double propPos;	// proportion positive
	
	// test that ut_r.decide produces a 80% of the time and b 20%:
	propPos = determineNTimes( N, ut_r, ESDecisionValue(), hd, dvMap.get( "myR", "a" ) );
	TS_ASSERT_DELTA( propPos, .8, LIM );
	
	// deterministic decision
	vector<string> vals;
	vals += "1","2";
	dvMap.add_decision_values( "i", vals );
	
	scnXml::HSESDecision ut_d_xml ("\
		i(1): a\
		i(2): b",
	    "d",	// decision
	    "i",	// depends
	    "a,b"	// values
	);
	ESDecisionRandom ut_d( dvMap, ut_d_xml );
	
	TS_ASSERT_EQUALS( ut_d.determine( dvMap.get( "i", "1" ), hd ), dvMap.get( "d", "a" ) );
	TS_ASSERT_EQUALS( ut_d.determine( ESDecisionValue(), hd ), ESDecisionValue() );	// void input and output
	TS_ASSERT_EQUALS( ut_d.determine( dvMap.get( "i", "2" ), hd ), dvMap.get( "d", "b" ) );
	
	scnXml::HSESDecision ut_void_xml ("\
		i(1): a\
		i(2): b\
		i(void): b",		// depending on a void input is illegal
	    "e",	// decision
	    "i",	// depends
	    "a,b"	// values
	);
	// Check we get an error (this isn't the best message, but it does what we need):
	TS_ASSERT_THROWS_EQUALS( ESDecisionRandom ut_e( dvMap, ut_void_xml ), const std::runtime_error &e, string(e.what()), "i(void) encountered: void is not an outcome of i" );
    }
    
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
	const int N = 10000;
	const double LIM = .02;
	double propPos;	// proportion positive
	UnittestUtil::setTotalParasiteDensity( *whm, 0. );	// no parasites (so we test specificity)
	propPos = determineNTimes( N, d, dvMap.get( "test", "microscopy" ), hd, dvMap.get( "result", "negative" ) );
	TS_ASSERT_DELTA ( propPos, .75, LIM );
	propPos = determineNTimes( N, d, dvMap.get( "test", "RDT" ), hd, dvMap.get( "result", "negative" ) );
	TS_ASSERT_DELTA ( propPos, .942, LIM );
	TS_ASSERT_EQUALS( d.determine( dvMap.get( "test", "none" ), hd ), ESDecisionValue() );
	
	UnittestUtil::setTotalParasiteDensity( *whm, 80. );	// a few parasites (so we test sensitivity with 0-100 parasites)
	propPos = determineNTimes( N, d, dvMap.get( "test", "microscopy" ), hd, dvMap.get( "result", "positive" ) );
	TS_ASSERT_DELTA ( propPos, .75, LIM );
	propPos = determineNTimes( N, d, dvMap.get( "test", "RDT" ), hd, dvMap.get( "result", "positive" ) );
	TS_ASSERT_DELTA ( propPos, .539, LIM );
	TS_ASSERT_EQUALS( d.determine( dvMap.get( "test", "none" ), hd ), ESDecisionValue() );
    }
    
private:
    ESDecisionValueMap dvMap;
    WithinHostModel* whm;
};

#endif
