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
    ESDecisionTreeSuite () : whm( WithinHostModel::createWithinHostModel() ), hd( numeric_limits< double >::quiet_NaN(), *whm, NONE ) {}
    ~ESDecisionTreeSuite () {
	delete whm;
    }
    
    void setUp () {
	UnittestUtil::EmpiricalWHM_setup();
	// could seed random-number-generator, but shouldn't affect outcomes
	UnittestUtil::PkPdSuiteSetup (PkPd::PkPdModel::LSTM_PKPD);
	
	dvMap = new ESDecisionValueMap;
	hd.ageYears = numeric_limits< double >::quiet_NaN();
	hd.pgState = NONE;
    }
    void tearDown () {
	delete dvMap;
	UnittestUtil::PkPdSuiteTearDown ();
    }
    
    void testESDecisionValue () {
	vector<string> vals;
	vals += "1", "2", "3";
	dvMap->add_decision_values( "a", vals );
	
	// Assert getting unknown values/decisions fails:
	TS_ASSERT_THROWS_EQUALS( dvMap->get( "b", "1" ), const std::runtime_error &e, string(e.what()), "ESDecisionValueMap::get(): no decision b" );
	TS_ASSERT_THROWS_EQUALS( dvMap->get( "a", "4" ), const std::runtime_error &e, string(e.what()), "ESDecisionValueMap::get(): no value a(4)" );
	
	// check we can re-add the same decision, even when values are reordered:
	vals[0] = "3";
	vals[2] = "1";
	TS_ASSERT_THROWS_NOTHING( dvMap->add_decision_values( "a", vals ) );
	// but not if it's different:
	vals[1] = "5";
	TS_ASSERT_THROWS_EQUALS( dvMap->add_decision_values( "a", vals ), const std::runtime_error &e, string(e.what()), "CaseManagement: decision a's values don't match; expected value: 2" );
	vals += "2", "6";	// now "3", "5", "1", "2", "6"
	TS_ASSERT_THROWS_EQUALS( dvMap->add_decision_values( "a", vals ), const std::runtime_error &e, string(e.what()), "CaseManagement: decision a's values don't match; unexpected values: 5 6" );
	
	// Check a value can't occur twice:
	vals.clear();
	vals.resize(2,"9");	// "9","9"
	TS_ASSERT_THROWS_EQUALS( dvMap->add_decision_values( "c", vals ), const std::runtime_error &e, string(e.what()), "CaseManagement: decision c's value 9 in value list twice!" );
	
	vals[1] = "A";	// "9", "A"
	dvMap->add_decision_values( "c", vals );
	
	// check what outcomes we get:
	typedef ESDecisionValue::id_type id_type;
	// a should have values 0, 1, 2
	TS_ASSERT_EQUALS( dvMap->get( "a", "1" ).id, static_cast<id_type>(0) );
	TS_ASSERT_EQUALS( dvMap->get( "a", "2" ).id, static_cast<id_type>(1) );
	TS_ASSERT_EQUALS( dvMap->get( "a", "3" ).id, static_cast<id_type>(2) );
	// c should have values 0, 4
	TS_ASSERT_EQUALS( dvMap->get( "c", "9" ).id, static_cast<id_type>(0) );
	TS_ASSERT_EQUALS( dvMap->get( "c", "A" ).id, static_cast<id_type>(4) );
	
	ostringstream val_string;
	val_string << dvMap->format( dvMap->get( "a", "1" ) | dvMap->get( "c", "9" ) );
	TS_ASSERT_EQUALS( val_string.str(), "a(1), c(9)" );
    }
    
    /* Runs d.determine( input, hd ) N times
     * returns the proportion of these runs where the output equalled expectedOutput
     */
    double determineNTimes (int N, const ESDecisionTree& d, const ESDecisionValue input, const ESHostData& hd, const ESDecisionValue expectedOutput) {
	int nExpected = 0;
	for (int i = 0; i < N; ++i) {
	    if( d.determine( input & d.mask, hd ) == expectedOutput )
		++nExpected;
	}
	return double(nExpected) / double(N);
    }
    
    void testRandomP () {
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
	ESDecisionRandom ut_r( *dvMap, ut_r_xml );
	
	const int N = 10000;
	const double LIM = .02;
	double propPos;	// proportion positive
	
	// test that ut_r.decide produces a 80% of the time and b 20%:
	propPos = determineNTimes( N, ut_r, ESDecisionValue(), hd, dvMap->get( "myR", "a" ) );
	TS_ASSERT_DELTA( propPos, .8, LIM );
	
    }
    
    void testRandomDeterministic () {
	// deterministic decision
	vector<string> vals;
	vals += "1","2";
	dvMap->add_decision_values( "i", vals );
	
	scnXml::HSESDecision ut_d_xml ("\
		i(1)	:a\
		i ( 2 )	{ b }",
	    "ut_d",	// decision
	    "i",	// depends
	    "b,a"	// values
	);
	ESDecisionRandom ut_d( *dvMap, ut_d_xml );
	
	TS_ASSERT_EQUALS( ut_d.determine( dvMap->get( "i", "1" ) & ut_d.mask, hd ), dvMap->get( "ut_d", "a" ) );
	TS_ASSERT_EQUALS( ut_d.determine( dvMap->get( "i", "2" ) & ut_d.mask, hd ), dvMap->get( "ut_d", "b" ) );
    }
    
    void testRandomErrors () {
	vector<string> vals;
	vals += "1","2";
	dvMap->add_decision_values( "i", vals );
	dvMap->add_decision_values( "j", vals );
	
	scnXml::HSESDecision ut_bad_decis_xml ("\
		i(1): a\
		j(2): b",	// decision changes here (illegal & non-sensical)
	    "bad_decis",	// decision
	    "i",	// depends
	    "a,b"	// values
	);
	TS_ASSERT_THROWS_EQUALS(
	    ESDecisionRandom ut_e( *dvMap, ut_bad_decis_xml ),
	    const std::runtime_error &e,
	    string(e.what()),
	    "failed to parse tree for bad_decis; expecting: \"i\" here: \"j(2): b\""
	);
	
	scnXml::HSESDecision ut_unknown_input_value_xml ("\
		i(1): a\
		i(2): b\
		i(3): b",		// depending on an unknown input value is illegal
	    "unknown_input_value",	// decision
	    "i",	// depends
	    "a,b"	// values
	);
	TS_ASSERT_THROWS_EQUALS(
	ESDecisionRandom ut_e( *dvMap, ut_unknown_input_value_xml ),
	    const std::runtime_error &e,
	    string(e.what()),
	    "decision tree unknown_input_value: i(3) encountered: 3 is not an outcome of i"
	);
	
	scnXml::HSESDecision ut_undeclared_output_xml ("\
		i(1): a\
		i(2): b",		// output not declared (illegal)
	    "undeclared_output",	// decision
	    "i",	// depends
	    "a"	// values
	);
	TS_ASSERT_THROWS_EQUALS(
	ESDecisionRandom ut_e( *dvMap, ut_undeclared_output_xml ),
	    const std::runtime_error &e,
	    string(e.what()),
	    "ESDecisionValueMap::get(): no value undeclared_output(b)"
	);
	
	scnXml::HSESDecision ut_nondepends_xml ("\
		i(1) {\
		    j(1): a\
		    j(2): b\
		}\
		i(2) {\
		    j(1): b\
		    j(2): a\
		}\
	    ",
	    "nondepends",	// decision
	    "i",	// error: dependency j not listed
	    "a,b"	// values
	);
	TS_ASSERT_THROWS_EQUALS(
	    ESDecisionRandom ut_f( *dvMap, ut_nondepends_xml ),
	    const std::runtime_error &e,
	    string(e.what()),
	    "decision tree nondepends: j not listed as a dependency"
	);
    }
    
    void testUC2Test () {
	ESDecisionUC2Test d( *dvMap );
	hd.pgState = STATE_MALARIA;
	TS_ASSERT_EQUALS( d.determine( ESDecisionValue() & d.mask, hd ), dvMap->get( "case", "UC1" ) );
	hd.pgState = static_cast<State>( STATE_MALARIA | SECOND_CASE );
	TS_ASSERT_EQUALS( d.determine( ESDecisionValue() & d.mask, hd ), dvMap->get( "case", "UC2" ) );
    }
    
    void testAge5Test () {
	ESDecisionAge5Test d( *dvMap );
	hd.ageYears = 4.99;
	TS_ASSERT_EQUALS( d.determine( ESDecisionValue() & d.mask, hd ), dvMap->get( "age5Test", "under5" ) );
	hd.ageYears = 5.0;
	TS_ASSERT_EQUALS( d.determine( ESDecisionValue() & d.mask, hd ), dvMap->get( "age5Test", "over5" ) );
    }
    
    void testParasiteTest () {
	ESDecisionParasiteTest d( *dvMap );
	hd.pgState = STATE_MALARIA;
	const int N = 10000;
	const double LIM = .02;
	double propPos;	// proportion positive
	
	UnittestUtil::setTotalParasiteDensity( *whm, 0. );	// no parasites (so we test specificity)
	propPos = determineNTimes( N, d, dvMap->get( "test", "microscopy" ), hd, dvMap->get( "result", "negative" ) );
	TS_ASSERT_DELTA ( propPos, .75, LIM );
	propPos = determineNTimes( N, d, dvMap->get( "test", "RDT" ), hd, dvMap->get( "result", "negative" ) );
	TS_ASSERT_DELTA ( propPos, .942, LIM );
	TS_ASSERT_EQUALS( d.determine( dvMap->get( "test", "none" ) & d.mask, hd ), dvMap->get( "result", "none" ) );
	
	UnittestUtil::setTotalParasiteDensity( *whm, 80. );	// a few parasites (so we test sensitivity with 0-100 parasites)
	propPos = determineNTimes( N, d, dvMap->get( "test", "microscopy" ), hd, dvMap->get( "result", "positive" ) );
	TS_ASSERT_DELTA ( propPos, .75, LIM );
	propPos = determineNTimes( N, d, dvMap->get( "test", "RDT" ), hd, dvMap->get( "result", "positive" ) );
	TS_ASSERT_DELTA ( propPos, .539, LIM );
	TS_ASSERT_EQUALS( d.determine( dvMap->get( "test", "none" ) & d.mask, hd ), dvMap->get( "result", "none" ) );
    }
    
    void testESDecisionMap () {
	// Tests using multiple decisions and ESDecisionMap::determine()
	// We basially test all tree-execution behaviour at once here.
	
	xsd::cxx::tree::sequence<scnXml::HSESDecision, false> decisionSeq;
	decisionSeq.push_back(
	    scnXml::HSESDecision( "\
		p(.1): none\
		p(.9){\
		    age5Test(under5){p(.8):RDT p(.2):microscopy}\
		    age5Test(over5){p(.5):RDT p(.5):microscopy}\
		}",	// tree
		"test",	// decision
		"age5Test",	// depends
		"none,RDT,microscopy"	// values
	    )
	);
	decisionSeq.push_back(
	    scnXml::HSESDecision( "\
		test(none){ case(UC2):second case(UC1):normal }\
		test(microscopy){\
		    result(positive){ case(UC2):second case(UC1):normal }\
		    result(negative): minor\
		    result(none): error\
		}\
		test(RDT){\
		    result(positive){ case(UC2):second case(UC1):normal }\
		    result(negative): minor\
		    result(none): error\
		}\
		",	// tree
		"treatment",	// decision
		"test,result,case",	// depends
		"minor,normal,second,error"	// values
	    )
	);
	scnXml::HSESDecisions decisions;
	decisions.setDecision( decisionSeq );
	
	scnXml::HSESTreatments treatments;	// empty treatment list
	
	// Final CaseManagement element
	::scnXml::HSESCaseManagement xmlCM( decisions, treatments );
	
	// use uncomplicated tree with its extra tests
	ESDecisionMap dMap;
	dMap.initialize( xmlCM, false );
	
	hd.ageYears = 2;
	hd.pgState = static_cast<State>( STATE_MALARIA | SECOND_CASE );
	UnittestUtil::setTotalParasiteDensity( *whm, 4000. );	// lots of parasites
	
	const int N = 100000;	// number of times to sample
	const double LIM = .002;	// answer expected to be accurate to this limit
	// Note: LIM=.002 is on the verge of what worked; it may need to be increased.
	
	int nMinor = 0, nNormal = 0, nSecond = 0;
	ESDecisionValue mask = dMap.dvMap.getDecision( "treatment" ).first;
	ESDecisionValue minor = dMap.dvMap.get( "treatment", "minor" );
	ESDecisionValue normal = dMap.dvMap.get( "treatment", "normal" );
	ESDecisionValue second = dMap.dvMap.get( "treatment", "second" );
	
	for (int i = 0; i < N; ++i) {
	    ESDecisionValue outcome = mask & dMap.determine( hd );
	    if( outcome == minor ) nMinor++;
	    else if( outcome == normal ) nNormal++;
	    else if( outcome == second ) nSecond++;
	    else TS_FAIL( "unexpected treatment output (error?)" );
	}
	
	// route: tested & (RDT | microscopy) & positive
	TS_ASSERT_DELTA( nMinor / double(N), 0.9 * (0.8*0.008 + 0.2*0.1), LIM );
	// impossible
	TS_ASSERT_EQUALS( nNormal, 0 );
	// route: not tested | (tested & (RDT | microscopy) & positive)
	TS_ASSERT_DELTA( nSecond / double(N), 0.1 + 0.9 * (0.8*0.992 + 0.2*0.9), LIM );
    }
    
private:
    ESDecisionValueMap* dvMap;
    WithinHostModel* whm;
    ESHostData hd;
};

#endif
