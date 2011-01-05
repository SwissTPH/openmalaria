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
// Unittest for the LSTM drug model

#ifndef Hmod_LSTMPkPdSuite
#define Hmod_LSTMPkPdSuite

#include <cxxtest/TestSuite.h>
#include "PkPd/LSTMPkPdModel.h"
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include <limits>

using namespace OM;
using namespace OM::PkPd;

const double NaN = numeric_limits<double>::quiet_NaN();

class LSTMPkPdSuite : public CxxTest::TestSuite
{
public:
    void setUp () {
	UnittestUtil::PkPdSuiteSetup(PkPdModel::LSTM_PKPD);
	proxy = new LSTMPkPdModel ();
	proteome_ID = 0;		// 0 should work; we definitely don't want random allocation
    }
    void tearDown () {
	delete proxy;
	UnittestUtil::PkPdSuiteTearDown();
    }
    
    void testNone () {
	TS_ASSERT_EQUALS (proxy->getDrugFactor (proteome_ID), 1.0);
    }
    
    void testOral () {
	proxy->medicate ("MF", 3000, 0, NaN, 21);
	TS_ASSERT_APPROX (proxy->getDrugFactor (proteome_ID), 0.03564073617400945);
    }
    
    void testOralHalves () {	// the point being: check it can handle two doses at the same time-point correctly
        //Note: normally NaN is used for duration, but 0 should give same result
	proxy->medicate ("MF", 1500, 0, 0, 21);
	proxy->medicate ("MF", 1500, 0, 0, 21);
	TS_ASSERT_APPROX (proxy->getDrugFactor (proteome_ID), 0.03564073617400945);
    }
    
    void testOralSplit () {
	proxy->medicate ("MF", 3000, 0, NaN, 21);
	proxy->medicate ("MF", 0, 0.5, NaN, 21);	// insert a second dose half way through the day: forces drug calculation to be split into half-days but shouldn't affect result
	TS_ASSERT_APPROX (proxy->getDrugFactor (proteome_ID), 0.03564073617400945);
    }
    
    void testOralDecayed () {
	proxy->medicate ("MF", 3000, 0, NaN, 21);
	proxy->decayDrugs ();
	TS_ASSERT_APPROX (proxy->getDrugFactor (proteome_ID), 0.03601694155274731);
    }
    
    void testOral2Doses () {
	proxy->medicate ("MF", 3000, 0, NaN, 21);
	proxy->decayDrugs ();
	proxy->medicate ("MF", 3000, 0, NaN, 21);
	TS_ASSERT_APPROX (proxy->getDrugFactor (proteome_ID), 0.03245158219000328);
    }
    
    // IV tests. MF may not be used as an IV drug, but we can still use it to test.
    void testIV () {
        // IV over whole day
        proxy->medicate ("MF", 50, 0, 1, 21);
        TS_ASSERT_APPROX (proxy->getDrugFactor(proteome_ID), 0.02212236680144665);
    }
    
    void testIVSplit (){
        // As above, but split into two doses
        proxy->medicate ("MF", 50, 0, 0.5, 21);
        proxy->medicate ("MF", 50, 0.5, 0.5, 21);
        TS_ASSERT_APPROX (proxy->getDrugFactor(proteome_ID), 0.02212236680144665);
    }
    
    void testCombined (){
        proxy->medicate ("MF", 50, 0, 0.5, 21);
        proxy->medicate ("MF", 1500, 0.5, NaN, 21);
        TS_ASSERT_APPROX (proxy->getDrugFactor(proteome_ID), 0.02714870841762252);
    }
    
    void testSimultaneous (){
        proxy->medicate ("MF", 1500, 0, NaN, 21);
        proxy->medicate ("MF", 50, 0, 0.5, 21);
        TS_ASSERT_APPROX (proxy->getDrugFactor(proteome_ID), 0.02715899967167571);
    }
    void testSimultaneousReversed (){
        // Note: IV dose registered first. Drug code must rearrange these to work correctly.
        proxy->medicate ("MF", 50, 0, 0.5, 21);
        proxy->medicate ("MF", 1500, 0, NaN, 21);
        TS_ASSERT_APPROX (proxy->getDrugFactor(proteome_ID), 0.02715899967167571);
    }
    
    LSTMPkPdModel *proxy;
    uint32_t proteome_ID;
};

#endif
