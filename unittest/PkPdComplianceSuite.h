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
        MQ_index = LSTMDrugType::findDrug( "MQ" );
    }
    void tearDown () {
        delete proxy;
        LSTMDrugType::clear();
    }

    void testMQ () {
        const double dose = 8.3 * bodymass;   // 8.3 mg/kg * 50 kg
        std::multimap<size_t, double> doses_for_day;
        doses_for_day.insert(make_pair(0, 0));
        doses_for_day.insert(make_pair(1, 0));
        doses_for_day.insert(make_pair(2, 0));
        const double drug_conc[] = { 0, 0.378440101, 0.737345129, 1.077723484,
                1.022091411, 0.969331065 };
        const double drug_factors[] = { 1, 0.031745814, 0.001007791, 0.000032,
                0.00000102, 3.22E-008 };
        typedef std::multimap<size_t,double>::const_iterator iter;
        for( size_t i = 0; i < 6; i++){
            // before update (after last step):
            //TODO: we think this should be 1 the first TWO steps if order is as we understand it:
            TS_ASSERT_APPROX_TOL (proxy->getDrugFactor (genotype), drug_factors[i], 1e-4, 1e-4);
//             double fac = proxy->getDrugFactor(genotype);
            // update (two parts):
            UnittestUtil::incrTime(sim::oneDay());
            proxy->decayDrugs();
            // after update:
            TS_ASSERT_APPROX_TOL (proxy->getDrugConc(MQ_index), drug_conc[i], 1e-4, 1e-4);
//             cout << "Conc, factor: " << proxy->getDrugConc(MQ_index) << ", " << fac << endl;
            // medicate (take effect on next update):
            pair<iter, iter> doses_tmp = doses_for_day.equal_range(i);
            for( iter it = doses_tmp.first; it != doses_tmp.second; it++){
                UnittestUtil::medicate( *proxy, MQ_index, dose, it->second,
                        numeric_limits< double >::quiet_NaN(), bodymass );
            }
        }
    }
    
private:
    LSTMModel *proxy;
    uint32_t genotype;
    double bodymass;
    size_t MQ_index;
};

#endif
