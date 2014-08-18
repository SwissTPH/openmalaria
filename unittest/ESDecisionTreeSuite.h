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
    
    double testMgPrescribed( scnXml::DecisionTree& dt, double age ){
        whm->age = age;
        UnittestUtil::clearMedicateQueue( whm->pkpd );
        TS_ASSERT_EQUALS( propTreatmentsNReps( 1, dt, *hd ), 1 );
        return UnittestUtil::getPrescribedMg( whm->pkpd );
    }
    
    void testDosing(){
        scnXml::DTTreatPKPD treat1( "sched1", "dosage1" );
        scnXml::DecisionTree simpleTreat;
        simpleTreat.getTreatPKPD().push_back( treat1 );
        
        // Test our dosing table. Set with a multiplier of 1 below 5 and 5 from 5.
        TS_ASSERT_DELTA( testMgPrescribed( simpleTreat, 0 ), 6, 1e-8 );
        TS_ASSERT_DELTA( testMgPrescribed( simpleTreat, 4.9 ), 6, 1e-8 );
        TS_ASSERT_DELTA( testMgPrescribed( simpleTreat, 5 ), 30, 1e-8 );
        TS_ASSERT_DELTA( testMgPrescribed( simpleTreat, 99 ), 30, 1e-8 );
    }
    
private:
    auto_ptr<WHMock> whm;
    auto_ptr<CMHostData> hd;
};

#endif
