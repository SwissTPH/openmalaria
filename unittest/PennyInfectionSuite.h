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

#ifndef Hmod_PennyInfectionSuite
#define Hmod_PennyInfectionSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"
#include "WithinHost/Infection/PennyInfection.h"
#include "util/random.h"
#include <limits>
#include <fstream>
#include <iomanip>

using namespace OM::WithinHost;

class PennyInfectionSuite : public CxxTest::TestSuite
{
public:
    void setUp () {
        m_rng.seed(1095);
        UnittestUtil::initTime(1);
        UnittestUtil::Infection_init_latentP_and_NaN ();
        PennyInfection::init();
        util::global_RNG.seed( 1095 );
        infection = new PennyInfection (m_rng, 0xFFFFFFFF);    // pkpdID (value) isn't important since we're not using drug model here
    }
    void tearDown () {
        delete infection;
    }
    
    void testThresholds(){
//         cout << setprecision(8) << infection->threshold_N
//             << '\t' << infection->threshold_C
//             << '\t' << infection->threshold_V << endl;
        TS_ASSERT_APPROX( infection->threshold_N, 8181.5227 );
        TS_ASSERT_APPROX( infection->threshold_C, 413.22176 );
        TS_ASSERT_APPROX( infection->threshold_V, 774.69253 );
    }
    
    static void readVector(std::vector<double>& vec, const char* file){
        ifstream istr(file);
        ETS_ASSERT( istr.is_open() );
        double val;
        while( istr >> val ){
            vec.push_back( val );
        }
        ETS_ASSERT( istr.eof() );
    }
    
    void testDensities(){
        vector<double> cirDens, seqDens;
        readVector(cirDens,"PennyCirDens.txt");
        readVector(seqDens,"PennySeqDens.txt");
//         ofstream outCir("PennyCirDens.txt");
//         ofstream outSeq("PennySeqDens.txt");
//         outCir << setprecision(8);
//         outSeq << setprecision(8);
        TS_ASSERT_EQUALS( cirDens.size(), seqDens.size() );
        
        bool extinct = false;
        int iterations=0;
        SimTime now = sim::ts0();
        do{
            extinct = infection->update(m_rng, 1.0, now, numeric_limits<double>::quiet_NaN());
            int ageDays = (now - infection->m_startDate - infection->s_latentP).inDays();
            while( ageDays < 0 ) ageDays += infection->delta_V; // special case encountered by unit test
            ETS_ASSERT_LESS_THAN( iterations, cirDens.size() );
            TS_ASSERT_APPROX( infection->getDensity(), cirDens[iterations] );
            TS_ASSERT_APPROX( infection->seqDensity(ageDays), seqDens[iterations] );
//             outCir << infection->getDensity() << endl;
//             outSeq << infection->seqDensity(ageDays) << endl;
            now += SimTime::oneDay();
            iterations+=1;
        }while(!extinct);
        TS_ASSERT_EQUALS( iterations, cirDens.size() );
    }
    
private:
    PennyInfection* infection;
    LocalRng m_rng;
};

#endif
