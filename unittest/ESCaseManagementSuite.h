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
#include <list>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'

using namespace OM::Clinical;
using namespace boost::assign;

class ESCaseManagementSuite : public CxxTest::TestSuite
{
public:
    //TODO: test tree execution (that all necessary decisions are evaluated and outputs conglomerated as expected)
    
    //TODO: test trees handle "void" output correctly
    
    void setUp () {
	vector<string> vals;
	vals += "extra","poor";
	dMap.dvMap.add_decision_values( "modQty", vals );
	vals.clear();
	vals += "0","5";
	dMap.dvMap.add_decision_values( "modD1", vals );
	vals.clear();
	vals += "B2";
	dMap.dvMap.add_decision_values( "modD2", vals );
	vals.clear();
	vals += "all", "selective";
	dMap.dvMap.add_decision_values( "modSTR", vals );
	
	// We need to add a "treatment" decision to make the thing work, though we won't be using it:
	xsd::cxx::tree::sequence<scnXml::HSESDecision, false> decisionSeq;
	decisionSeq.push_back(
	    scnXml::HSESDecision(
		"void",	// tree
		"treatment",	// decision
		"",	// depends
		"treatment1"	// values
	    )
	);
	scnXml::HSESDecisions decisions;
	decisions.setDecision( decisionSeq );
	
	// Base schedule
	xsd::cxx::tree::sequence<scnXml::Medicate, false> medicateSeq;
	medicateSeq.push_back( scnXml::Medicate( "A", 1000.0, 0 ) );
	medicateSeq.push_back( scnXml::Medicate( "B", 3000.0, 0 ) );
	medicateSeq.push_back( scnXml::Medicate( "B", 3000.0, 12 ) );
	scnXml::HSESTreatmentSchedule treatSched;
	treatSched.setMedicate( medicateSeq );
	
	// Modifiers
	xsd::cxx::tree::sequence<scnXml::HSESTreatmentModifierEffect, false> modQtySeq;
	modQtySeq.push_back( scnXml::HSESTreatmentModifierEffect( "extra", "A(2),B(1.3)" ) );
	modQtySeq.push_back( scnXml::HSESTreatmentModifierEffect( "poor", " A(0.5) , B( 0.2 ) " ) );
	scnXml::HSESTreatmentModifier modQty( "modQty" );
	modQty.setMultiplyQty( modQtySeq );
	
	xsd::cxx::tree::sequence<scnXml::HSESTreatmentModifierEffect, false> modD1Seq;
	modD1Seq.push_back( scnXml::HSESTreatmentModifierEffect( "0", "A(0),B(0)" ) );
	modD1Seq.push_back( scnXml::HSESTreatmentModifierEffect( "5", "A(5),B(5)" ) );
	scnXml::HSESTreatmentModifier modD1( "modD1" );
	modD1.setDelay( modD1Seq );
	
	xsd::cxx::tree::sequence<scnXml::HSESTreatmentModifierEffect, false> modD2Seq;
	modD2Seq.push_back( scnXml::HSESTreatmentModifierEffect( "B2", "B(2 ),A(0) " ) );	// note: backwards
	scnXml::HSESTreatmentModifier modD2( "modD2" );
	modD2.setDelay( modD2Seq );
	
	xsd::cxx::tree::sequence<scnXml::HSESTreatmentModifierEffect, false> modSTRSeq;
	modSTRSeq.push_back( scnXml::HSESTreatmentModifierEffect( "all", "A(0-100),B( 0-100)" ) );
	modSTRSeq.push_back( scnXml::HSESTreatmentModifierEffect( "selective", "A(0-1000 ),B(2-100)" ) );
	scnXml::HSESTreatmentModifier modSTR( "modSTR" );
	modSTR.setSelectTimeRange( modSTRSeq );
	
	xsd::cxx::tree::sequence<scnXml::HSESTreatmentModifier, false> modifierSeq;
	modifierSeq.push_back( modQty );
	modifierSeq.push_back( modD1 );
	modifierSeq.push_back( modD2 );
	modifierSeq.push_back( modSTR );
	
	// Treatment
	scnXml::HSESTreatment treatment1( treatSched, "treatment1" );
	treatment1.setModifier( modifierSeq );
	
	scnXml::HSESTreatments treatments;
	xsd::cxx::tree::sequence<scnXml::HSESTreatment, false> treatmentSeq( 1, treatment1 );
	treatments.setTreatment( treatmentSeq );
	
	// Final CaseManagement element
	::scnXml::HSESCaseManagement xmlCM( decisions, treatments );
	
	// use complicated tree, because it doesn't add so many unwanted decisions
	dMap.initialize( xmlCM, true );
    }
    
    void testTreatments() {
	// Note: use ETS_..., not TS_..., when assertion should throw if false (prevent dangerous operations)
	
	const ESTreatmentSchedule *sched;
	sched = dMap.getSchedule( ESDecisionValue() );	// void input should resolve an empty schedule
	ETS_ASSERT( sched != NULL );
	list<MedicateData> medQueue;
	sched->apply( medQueue );
	TS_ASSERT_EQUALS( medQueue.size(), 0u );
	
	ESDecisionValue treatment1 = dMap.dvMap.get( "treatment", "treatment1" );	// has 3 treatments; A at time 0 and B at times 0,12
	treatment1 |= dMap.dvMap.get( "modQty", "poor" );	// reduce quantities
	treatment1 |= dMap.dvMap.get( "modD1", "5" );	// delay by 5 hours
	treatment1 |= dMap.dvMap.get( "modD2", "B2" );	// delay B by 2 hours
	treatment1 |= dMap.dvMap.get( "modSTR", "selective" );	// A-0 and B-12 should be kept (B-0 shouldn't, since delays should be added after selection)
	
	sched = dMap.getSchedule( treatment1 );
	ETS_ASSERT( sched != NULL );
	
	medQueue.clear();
	sched->apply( medQueue );
	ETS_ASSERT_EQUALS( medQueue.size(), 2u );
	
	const MedicateData& md1 = medQueue.front();
	TS_ASSERT_EQUALS( md1.abbrev, "A" );
	TS_ASSERT_EQUALS( md1.qty, 500.0 );
	TS_ASSERT_DELTA( md1.time, 5.0/24.0, 1.0e-10 );
	medQueue.pop_front();
	const MedicateData& md2 = medQueue.front();
	TS_ASSERT_EQUALS( md2.abbrev, "B" );
	TS_ASSERT_EQUALS( md2.qty, 600.0 );
	TS_ASSERT_DELTA( md2.time, 19.0/24.0, 1.0e-10 );
    }
    
private:
    ESDecisionMap dMap;
};

#endif
