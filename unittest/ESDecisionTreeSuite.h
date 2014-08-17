/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
 
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
#include "util/random.h"
#include "UnittestUtil.h"
#include "WHMock.h"
#include <limits>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

using namespace OM::Clinical;
using namespace OM::WithinHost;
using namespace boost::assign; // bring 'operator+=()' into scope
using UnitTest::WHMock;

class ESDecisionTreeSuite : public CxxTest::TestSuite
{
public:
    ESDecisionTreeSuite () :
            whm(0), hd(0)
    {
        UnittestUtil::initSurveys();
    }
    ~ESDecisionTreeSuite () {}
    
    void setUp () {
	// Note: cannot create whm in constructor, since it uses random number
	// generator which is initialized after constructor runs.
	util::random::seed (83);	// seed is unimportant, but must be fixed
	UnittestUtil::EmpiricalWHM_setup();     // use a 1-day-TS model
        whm.reset( new WHMock() );
        ETS_ASSERT( whm.get() != 0 );
	hd.reset( new CMHostData( numeric_limits< double >::quiet_NaN(), *whm.get(), Episode::NONE ) );

	UnittestUtil::EmpiricalWHM_setup();
	// could seed random-number-generator, but shouldn't affect outcomes
	UnittestUtil::PkPdSuiteSetup (PkPd::PkPdModel::LSTM_PKPD);
	
	hd->ageYears = numeric_limits< double >::quiet_NaN();
	hd->pgState = Episode::NONE;
    }
    void tearDown () {
        PkPd::LSTMDrugType::clear();
        PkPd::LSTMTreatments::clear();
    }
    
    /* Runs the decision tree N times
     * returns the proportion of these runs where the output was any treatment
     */
    double propTreatmentsNReps (int N, const scnXml::DecisionTree& dt, const CMHostData& hd) {
        auto_ptr<CMDecisionTree> cmdt = CMDecisionTree::create( dt );
        
	whm->nTreatments = 0;
	for (int i = 0; i < N; ++i) {
	    cmdt->exec( hd );
	}
	return double(whm->nTreatments) / double(N);
    }
    
    void testRandomP () {
	CMHostData hd(
		numeric_limits< double >::quiet_NaN(),
		*whm,
		Episode::NONE
	);
	
	// random decision
        // option (a) is to treat, option (b) is to do nothing
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        
        scnXml::Outcome o1r2( 0.9 ), o2r2( 0.1 );
        o1r2.getTreatPKPD().push_back( treat1 );
        o2r2.setNoAction( scnXml::DTNoAction() );
        scnXml::DTRandom r2;
        r2.getOutcome().push_back( o1r2 );
        r2.getOutcome().push_back( o2r2 );
        
        scnXml::Outcome o1r3( 0.7 ), o2r3( 0.3 );
        o1r3.getTreatPKPD().push_back( treat1 );
        o2r3.setNoAction( scnXml::DTNoAction() );
        scnXml::DTRandom r3;
        r3.getOutcome().push_back( o1r3 );
        r3.getOutcome().push_back( o2r3 );
        
        scnXml::Outcome o1r1( 0.5 ), o2r1( 0.5 );
        o1r1.setRandom( r2 );
        o2r1.setRandom( r3 );
        
        scnXml::DTRandom r1;
        r1.getOutcome().push_back( o1r1 );
        r1.getOutcome().push_back( o2r1 );
        
        scnXml::DecisionTree dt;
        dt.setRandom( r1 );
	
	const int N = 10000;
	const double LIM = .02;
	
        // test that dt.exec chooses to treat 80% and no action 20% of the time:
	TS_ASSERT_DELTA( propTreatmentsNReps( N, dt, hd ), 0.8, LIM );
    }
    
    //TODO: test dosing is correct
    
    void testUC2Test () {
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back( treat1 );
        scnXml::DecisionTree noAction;
        noAction.setNoAction( scnXml::DTNoAction() );
        
        scnXml::DTCaseType ct( simpleTreat,  // first line: simple treatment
                            noAction );     // second line: no action
        scnXml::DecisionTree dt;
        dt.setCaseType( ct );
        
	hd->pgState = static_cast<Episode::State>( Pathogenesis::STATE_MALARIA );
	TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt, *hd ), 1 );
	hd->pgState = static_cast<Episode::State>( Pathogenesis::STATE_MALARIA |
                Episode::SECOND_CASE );
	TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt, *hd ), 0 );
    }
    
    void testParasiteTest () {
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back( treat1 );
        scnXml::DecisionTree noAction;
        noAction.setNoAction( scnXml::DTNoAction() );
        
        scnXml::DTDiagnostic microscopy( simpleTreat,  // positive: simple treatment
                            noAction,     // negative: no action
                            "microscopy" );     // type of diagnostic
        scnXml::DecisionTree dt_mic;
        dt_mic.setDiagnostic( microscopy );
        
        scnXml::DTDiagnostic rdt( simpleTreat,  // positive: simple treatment
                            noAction,     // negative: no action
                            "RDT" );     // type of diagnostic
        scnXml::DecisionTree dt_rdt;
        dt_rdt.setDiagnostic( rdt );
        
        hd->pgState = static_cast<Episode::State>( Pathogenesis::STATE_MALARIA );
	const int N = 20000;
	const double LIM = .02;
	
	whm->totalDensity = 0.0;	// no parasites (so we test specificity)
	TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_mic, *hd ), 1 - 0.75, LIM );
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_rdt, *hd ), 1 - 0.942, LIM );
	
        whm->totalDensity = 80.0;	// a few parasites
	TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_mic, *hd ), 0.85, LIM );
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_rdt, *hd ), 0.63769, LIM );
	
        whm->totalDensity = 2000.0;	// lots of parasites
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_mic, *hd ), 0.99257, LIM );
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_rdt, *hd ), 0.99702, LIM );
    }
    
#if 0
    void xtestESDecisionMap () {
	// Tests using multiple decisions and ESDecisionMap::determine()
	// We basially test all tree-execution behaviour at once here.
	
	xsd::cxx::tree::sequence<scnXml::HSESDecision, false> decisionSeq;
	decisionSeq.push_back(
	    scnXml::HSESDecision( "\
		age(0-5): under5\
		age(5-inf): over5",	// tree
		"age5",	// decision
		"age",	// depends
		"under5,over5"	// values
	    )
	);
	decisionSeq.push_back(
	    scnXml::HSESDecision( "\
		p(.1): none\
		p(.9){\
		    age5(under5){p(.8):RDT p(.2):microscopy}\
		    age5(over5){p(.5):RDT p(.5):microscopy}\
		}",	// tree
		"test",	// decision
		"age5,p",	// depends
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
	dMap.initialize( xmlCM, ESDecisionMap::Uncomplicated, false );
	
	hd->ageYears = 2;
	hd->pgState = static_cast<Episode::State>( Pathogenesis::STATE_MALARIA | Episode::SECOND_CASE );
	whm->totalDensity = 4000.0;	// lots of parasites
	
	const int N = 100000;	// number of times to sample
	const double LIM = .002;	// answer expected to be accurate to this limit
	// Note: LIM=.002 is on the verge of what worked; it may need to be increased.
	
	int nMinor = 0, nNormal = 0, nSecond = 0;
	ESDecisionValue mask = dMap.dvMap.getDecisionMask( "treatment" );
	ESDecisionValue minor = dMap.dvMap.get( "treatment", "minor" );
	ESDecisionValue normal = dMap.dvMap.get( "treatment", "normal" );
	ESDecisionValue second = dMap.dvMap.get( "treatment", "second" );
	
	for (int i = 0; i < N; ++i) {
	    ESDecisionValue outcome = mask & dMap.determine( *hd );
	    if( outcome == minor ) nMinor++;
	    else if( outcome == normal ) nNormal++;
	    else if( outcome == second ) nSecond++;
	    else TS_FAIL( "unexpected treatment output (error?)" );
	}
	
	// route: tested & (RDT | microscopy) & positive
	TS_ASSERT_DELTA( nMinor / double(N), 0.9 * (0.8*0.012 + 0.2*0.004), LIM );
	// impossible
	TS_ASSERT_EQUALS( nNormal, 0 );
	// route: not tested | (tested & (RDT | microscopy) & positive)
	TS_ASSERT_DELTA( nSecond / double(N), 0.1 + 0.9 * (0.8*0.988 + 0.2*0.996), LIM );
    }
#endif
    
private:
    auto_ptr<WHMock> whm;
    auto_ptr<CMHostData> hd;
};

#endif
