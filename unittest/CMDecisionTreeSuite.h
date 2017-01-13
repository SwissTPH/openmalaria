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

#ifndef Hmod_CMDecisionTreeSuite
#define Hmod_CMDecisionTreeSuite

#include <cxxtest/TestSuite.h>
#include "Clinical/CMDecisionTree.h"
#include "util/random.h"
#include "UnittestUtil.h"
#include "WHMock.h"
#include <limits>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

using namespace OM::Clinical;
using namespace OM::WithinHost;
using namespace boost::assign; // bring 'operator+=()' into scope
using UnitTest::WHMock;

class CMDecisionTreeSuite : public CxxTest::TestSuite
{
public:
    CMDecisionTreeSuite () :
            whm(0), hd(0) {}
    ~CMDecisionTreeSuite () {}
    
    void setUp () {
        // Note: cannot create whm in constructor, since it uses random number
        // generator which is initialized after constructor runs.
        util::random::seed (83);	// seed is unimportant, but must be fixed
        UnittestUtil::initTime(5);
        UnittestUtil::initSurveys();
        UnittestUtil::setDiagnostics();

        UnittestUtil::EmpiricalWHM_setup();
        whm.reset( new WHMock() );
        ETS_ASSERT( whm.get() != 0 );

        human.reset( UnittestUtil::createHuman(SimTime::zero()).release() );
        ETS_ASSERT( human.get() != 0 );
        UnittestUtil::setHumanWH( *human, whm.get() );

        hd.reset( new CMHostData( *human, 21 /* any real age will do for most tests */,
                                      Episode::NONE ) );

        // could seed random-number-generator, but shouldn't affect outcomes
        UnittestUtil::PkPdSuiteSetup ();
    }
    void tearDown () {
        human.reset();
        whm.reset();
        PkPd::LSTMDrugType::clear();
        PkPd::LSTMTreatments::clear();
    }
    
    /* Runs the decision tree N times
     * returns the proportion of these runs where the output was any treatment
     */
    double propTreatmentsNReps (int N, const scnXml::DecisionTree& dt) {
        const CMDecisionTree& cmdt = CMDecisionTree::create( dt, true );
        
	whm->nTreatments = 0;
        int secondCounter = 0;
	for (int i = 0; i < N; ++i) {
	    secondCounter += (cmdt.exec( *hd ).treated ? 1 : 0);
	}
	TS_ASSERT_EQUALS( whm->nTreatments, secondCounter );
	return double(whm->nTreatments) / double(N);
    }
    
    void testRandomP () {
	// random decision
        // option (a) is to treat, option (b) is to do nothing
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        
        scnXml::Outcome o1r2( 0.9 ), o2r2( 0.1 );
        o1r2.getTreatPKPD().push_back( treat1 );
        o2r2.setNoTreatment( scnXml::DTNoTreatment() );
        scnXml::DTRandom r2;
        r2.getOutcome().push_back( o1r2 );
        r2.getOutcome().push_back( o2r2 );
        
        scnXml::Outcome o1r3( 0.7 ), o2r3( 0.3 );
        o1r3.getTreatPKPD().push_back( treat1 );
        o2r3.setNoTreatment( scnXml::DTNoTreatment() );
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
	TS_ASSERT_DELTA( propTreatmentsNReps( N, dt ), 0.8, LIM );
    }
    
    void testUC2Test () {
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back( treat1 );
        scnXml::DecisionTree noAction;
        noAction.setNoTreatment( scnXml::DTNoTreatment() );
        
        scnXml::DTCaseType ct( simpleTreat,  // first line: simple treatment
                            noAction );     // second line: no action
        scnXml::DecisionTree dt;
        dt.setCaseType( ct );
        
	hd->pgState = static_cast<Episode::State>( Pathogenesis::STATE_MALARIA );
	TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );
	hd->pgState = static_cast<Episode::State>( Pathogenesis::STATE_MALARIA |
                Episode::SECOND_CASE );
	TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 0 );
    }
    
    void testParasiteTest () {
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back( treat1 );
        scnXml::DecisionTree noAction;
        noAction.setNoTreatment( scnXml::DTNoTreatment() );
        
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
	TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_mic ), 1 - 0.75, LIM );
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_rdt ), 1 - 0.942, LIM );
	
        whm->totalDensity = 80.0;	// a few parasites
	TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_mic ), 0.96815, LIM );
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_rdt ), 0.90599, LIM );
	
        whm->totalDensity = 2000.0;	// lots of parasites
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_mic ), 0.99867, LIM );
        TS_ASSERT_DELTA( propTreatmentsNReps( N, dt_rdt ), 0.99584, LIM );
    }
    
    void testAgeSwitch(){
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DTAge ageSwitch;
        
        scnXml::Age treatYoung( 0.0 );  // ages 0 to 2.5
        treatYoung.getTreatPKPD().push_back( treat1 );
        ageSwitch.getAge().push_back( treatYoung );
        
        scnXml::Age noTreat( 2.5 );     // ages 2.5 to 50
        noTreat.setNoTreatment( scnXml::DTNoTreatment() );
        ageSwitch.getAge().push_back( noTreat );
        
        scnXml::Age treatOlder( 50.0 );  // ages 50+
        treatOlder.getTreatPKPD().push_back( treat1 );
        ageSwitch.getAge().push_back( treatOlder );
        
        scnXml::DecisionTree dt;
        dt.setAge( ageSwitch );
        
        hd->ageYears = 1;
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );
        hd->ageYears = 2.5;
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 0 );
        hd->ageYears = 50;
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );
        hd->ageYears = 1e6;     // there is no upper bound, so this should work!
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );
    }
    
    void testSimpleTreat(){
        TS_ASSERT_EQUALS( whm->lastTimeLiver, SimTime::never() );
        TS_ASSERT_EQUALS( whm->lastTimeBlood, SimTime::never() );
        
        scnXml::DTTreatSimple treat1( "0t", "1t" );   // 0 time steps liver, 1 blood (using 1 day TS)
        scnXml::DecisionTree dt1;
        dt1.setTreatSimple( treat1 );
        
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt1 ), 1 );
        TS_ASSERT_EQUALS( whm->lastTimeLiver.inDays(), 0 );
        TS_ASSERT_EQUALS( whm->lastTimeBlood.inDays(), 1*5 );
        
        scnXml::DTTreatSimple treat2( "15d", "-1t" );
        scnXml::DecisionTree dt2;
        dt2.setTreatSimple( treat2 );
        
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt2 ), 1 );
        TS_ASSERT_EQUALS( whm->lastTimeLiver.inDays(), 3*5 );
        TS_ASSERT_EQUALS( whm->lastTimeBlood.inDays(), -1*5 );
    }
    
    double runAndGetMgPrescribed( scnXml::DecisionTree& dt, double age ){
        hd->ageYears = age;
        UnittestUtil::clearMedicateQueue( whm->pkpd );
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt ), 1 );
        return UnittestUtil::getPrescribedMg( whm->pkpd );
    }
    
    void testDosing(){
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DecisionTree dt1;
        dt1.getTreatPKPD().push_back( treat1 );
        
        // Test our dosing table. Set with a multiplier of 1 below 5 and 5 from 5.
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt1, 0 ), 6, 1e-8 );
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt1, 4.9 ), 6, 1e-8 );
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt1, 5 ), 30, 1e-8 );
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt1, 99 ), 30, 1e-8 );
        
        scnXml::DTTreatPKPD treat2( "sched2", "dosage1" );
        scnXml::DecisionTree dt2;
        dt2.getTreatPKPD().push_back( treat2 );
        
        // Test our dosing table. Set with a multiplier of 1 below 5 and 5 from 5.
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt2, 0 ), 7, 1e-8 );
        TS_ASSERT_DELTA( runAndGetMgPrescribed( dt2, 99 ), 35, 1e-8 );
    }
    
private:
    auto_ptr<Host::Human> human;
    auto_ptr<WHMock> whm;
    auto_ptr<CMHostData> hd;
};

#endif
