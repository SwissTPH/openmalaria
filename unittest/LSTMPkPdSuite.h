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
// Unittest for the LSTM drug model

#ifndef Hmod_LSTMPkPdSuite
#define Hmod_LSTMPkPdSuite

#include <cxxtest/TestSuite.h>
#include "PkPd/LSTMModel.h"
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include <limits>

using namespace OM;
using namespace OM::PkPd;

const double NaN = numeric_limits<double>::quiet_NaN();

class LSTMPkPdSuite : public CxxTest::TestSuite
{
public:
    LSTMPkPdSuite() :
            proxy(0)
    {
        genotype = 0;                // 0 should work; we definitely don't want random allocation
        // This is what a previously-used weight distribution gave us,
        // and is good enough for testing purposes:
        massAt21 = 55.4993;
    }
    
    void setUp () {
        UnittestUtil::initTime(1);
	UnittestUtil::PkPdSuiteSetup();
	proxy = new LSTMModel ();
        MQ_index = LSTMDrugType::findDrug( "MQ" );
    }
    void tearDown () {
	delete proxy;
        LSTMDrugType::clear();
    }
    
    void testNone () {
	TS_ASSERT_EQUALS (proxy->getDrugFactor (genotype), 1.0);
    }
    
    void testOral () {
	UnittestUtil::medicate( *proxy, MQ_index, 3000, 0, NaN, massAt21 );
	TS_ASSERT_APPROX (proxy->getDrugFactor (genotype), 0.03174563638523168);
    }
    
    void testOralHalves () {	// the point being: check it can handle two doses at the same time-point correctly
        //Note: normally NaN is used for duration, but 0 should give same result
	UnittestUtil::medicate( *proxy, MQ_index, 1500, 0, 0, massAt21 );
	UnittestUtil::medicate( *proxy, MQ_index, 1500, 0, 0, massAt21 );
	TS_ASSERT_APPROX (proxy->getDrugFactor (genotype), 0.03174563638523168);
    }
    
    void testOralSplit () {
	UnittestUtil::medicate( *proxy, MQ_index, 3000, 0, NaN, massAt21 );
	UnittestUtil::medicate( *proxy, MQ_index, 0, 0.5, NaN, massAt21 );	// insert a second dose half way through the day: forces drug calculation to be split into half-days but shouldn't affect result
	TS_ASSERT_APPROX (proxy->getDrugFactor (genotype), 0.03174563639140275);
    }
    
    void testOralDecayed () {
	UnittestUtil::medicate( *proxy, MQ_index, 3000, 0, NaN, massAt21 );
	proxy->decayDrugs ();
	TS_ASSERT_APPROX (proxy->getDrugFactor (genotype), 0.03174563639501896);
    }
    
    void testOral2Doses () {
	UnittestUtil::medicate( *proxy, MQ_index, 3000, 0, NaN, massAt21 );
	proxy->decayDrugs ();
	UnittestUtil::medicate( *proxy, MQ_index, 3000, 0, NaN, massAt21 );
	TS_ASSERT_APPROX (proxy->getDrugFactor (genotype), 0.03174563637686205);
    }
    
    // IV tests. MQ may not be used as an IV drug, but our code doesn't care.
    void testIVEquiv () {
        // As duration tends to zero, factor should tend to that for an oral
        // dose. Code uses duration!=0 to enable IV mode so we use a small value;
        // if this is too small, however, the gsl_integration_qag function complains
        // it cannot reach the requested tolerance.
        UnittestUtil::medicate( *proxy, MQ_index, 3000/massAt21, 0, 1e-6, massAt21 );
        TS_ASSERT_APPROX (proxy->getDrugFactor (genotype), 0.03174563760775995);
    }
    
    void testIV () {
        // IV over whole day
        UnittestUtil::medicate( *proxy, MQ_index, 50, 0, 1, massAt21 );
        TS_ASSERT_APPROX (proxy->getDrugFactor(genotype), 0.03308874286174752);
    }
    
    void testIVSplit (){
        // As above, but split into two doses
        UnittestUtil::medicate( *proxy, MQ_index, 25, 0, 0.5, massAt21 );
        UnittestUtil::medicate( *proxy, MQ_index, 25, 0.5, 0.5, massAt21 );
        TS_ASSERT_APPROX (proxy->getDrugFactor(genotype), 0.03308874286174755);
    }
    
    void testCombined (){
        UnittestUtil::medicate( *proxy, MQ_index, 50, 0, 0.5, massAt21 );
        UnittestUtil::medicate( *proxy, MQ_index, 1500, 0.5, NaN, massAt21 );
        TS_ASSERT_APPROX (proxy->getDrugFactor(genotype), 0.03241010934374807);
    }
    
    void testSimultaneous (){
        UnittestUtil::medicate( *proxy, MQ_index, 1500, 0, NaN, massAt21 );
        UnittestUtil::medicate( *proxy, MQ_index, 50, 0, 0.5, massAt21 );
        TS_ASSERT_APPROX (proxy->getDrugFactor(genotype), 0.03174563640637330);
    }
    void testSimultaneousReversed (){
        // Note: IV dose registered first. Drug code must rearrange these to work correctly.
        UnittestUtil::medicate( *proxy, MQ_index, 50, 0, 0.5, massAt21 );
        UnittestUtil::medicate( *proxy, MQ_index, 1500, 0, NaN, massAt21 );
        TS_ASSERT_APPROX (proxy->getDrugFactor(genotype), 0.03174563640637330);
    }
    
private:
    LSTMModel *proxy;
    uint32_t genotype;
    double massAt21;
    size_t MQ_index;
};

#endif
