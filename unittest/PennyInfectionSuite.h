/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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
        TimeStep::init( 1, 90.0 );
        UnittestUtil::Infection_init_1day ();
        PennyInfection::init();
        TimeStep::simulation = TimeStep(7);     // value shouldn't be important
        util::random::seed( 1095 );
        infection = new PennyInfection (0xFFFFFFFF);    // pkpdID (value) isn't important since we're not using drug model here
    }
    void tearDown () {
        delete infection;
        util::random::seed(0);  // make sure nothing else uses this seed/reports
    }
    
    void testThresholds(){
        TS_ASSERT_APPROX( infection->threshold_N, 1162659.874442960284000 );
        TS_ASSERT_APPROX( infection->threshold_C, 220.766169197537423 );
        TS_ASSERT_APPROX( infection->threshold_V, 778.491750025389886 );
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
        TS_ASSERT_EQUALS( cirDens.size(), seqDens.size() );
        
        bool extinct = false;
        int iterations=0;
        do{
            extinct = infection->update(1.0);
            ETS_ASSERT_LESS_THAN( iterations, cirDens.size() );
            TS_ASSERT_APPROX( infection->getDensity(), cirDens[iterations] );
            TS_ASSERT_APPROX( infection->seqDensity(), seqDens[iterations] );
            TimeStep::simulation += TimeStep(1);
            iterations+=1;
        }while(!extinct);
        TS_ASSERT_EQUALS( iterations, cirDens.size() );
    }
    
private:
    PennyInfection* infection;
};

#endif
