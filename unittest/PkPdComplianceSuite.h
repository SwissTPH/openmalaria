/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
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

#ifndef Hmod_PkPdComplianceSuite
#define Hmod_PkPdComplianceSuite

#include <cxxtest/TestSuite.h>
#include "PkPd/LSTMModel.h"
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include <limits>

using namespace OM;
using namespace OM::PkPd;

#   define TS_ASSERT_APPROX6(x,y) ___TS_ASSERT_APPROX(__FILE__,__LINE__,x,y, 1e-6, 1e-6)

class PkPdComplianceSuite : public CxxTest::TestSuite
{
public:
    PkPdComplianceSuite() :
            proxy(0)
    {
        genotype = 0;                // 0 should work; we definitely don't want random allocation
        bodymass = 50;
	//List of drugs
	//singleCompartmentDrugs
    }
    
    void setUp () {
        UnittestUtil::initTime(1);
	UnittestUtil::PkPdSuiteSetup();
	proxy = new LSTMModel ();
        MF_index = LSTMDrugType::findDrug( "MF" );
    }
    void tearDown () {
        delete proxy;
        LSTMDrugType::clear();
    }

    void testMF () {
        const double drug_factors[] = { 1, 0.031745814, 0.001007791, 0.000032,
            0.00000102, 3.22E-008 };
        TS_ASSERT_APPROX6 (proxy->getDrugFactor (genotype), drug_factors[0]);
        UnittestUtil::incrTime(sim::oneDay());
        proxy->decayDrugs();
        UnittestUtil::medicate( *proxy, MF_index, 8300, 0, NaN, bodymass );
        TS_ASSERT_APPROX6 (proxy->getDrugFactor (genotype), drug_factors[1]);
    }
    
private:
    LSTMModel *proxy;
    uint32_t genotype;
    double bodymass;
    size_t MF_index;
};

#endif
