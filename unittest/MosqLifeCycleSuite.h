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

#ifndef Hmod_MosqLifeCycleSuite
#define Hmod_MosqLifeCycleSuite

#include <cxxtest/TestSuite.h>
#include "UnittestUtil.h"
#include "ExtraAsserts.h"

#include "Transmission/Anopheles/MosqTransmission.h"
#include "Transmission/Anopheles/LCEmergence.h"
#include "Transmission/Anopheles/ResourceFitter.h"
#include "schema/entomology.h"
#include "util/errors.h"

#include <gsl/gsl_vector.h>
#include <fstream>
#include <limits>


namespace OM { namespace Transmission { namespace Anopheles {
    // this is a private function of the ResourceFitter module not declared in
    // the header — but we can still test it
    void vector_scale_length( const gsl_vector *source, gsl_vector *target );
} } }
using namespace OM::Transmission::Anopheles;

// just some different constants
const double MLCS_a=5.23e-5, MLCS_b=9.4e12, MLCS_c=9.32, MLCS_d=9.34243e-2, MLCS_e=141.23;
const size_t YEAR_LEN = 365;
const double P_A = 0.685785;
const double P_df = 0.195997;
const double initNvFromSv =47.619, initOvFromSv = 3.71429;
const double myP_dif[] = { 0.0208892, 0.0211708, 0.0211384, 0.0207101, 0.020627, 0.020583, 0.0204871, 0.0202843, 0.0202156, 0.0200973, 0.0199916, 0.0199837, 0.0199941, 0.0198548, 0.0197392, 0.0196689, 0.0196281, 0.0196362, 0.0195606, 0.0196321, 0.0196513, 0.0196813, 0.0197236, 0.019749, 0.0198011, 0.0199042, 0.0199767, 0.0200928, 0.0203577, 0.020535, 0.0206487, 0.0207777, 0.020979, 0.0211049, 0.0212751, 0.0213585, 0.0213085, 0.0215146, 0.0215828, 0.021681, 0.0217439, 0.0218348, 0.0218485, 0.0218434, 0.0219025, 0.0219072, 0.0218668, 0.0218536, 0.0218229, 0.0217871, 0.0217856, 0.0217375, 0.0216027, 0.0214194, 0.0211886, 0.02102, 0.0210537, 0.0210596, 0.020977, 0.0208403, 0.0209002, 0.0209218, 0.0208807, 0.020998, 0.0209654, 0.0209916, 0.0210443, 0.0211605, 0.0211625, 0.0211092, 0.0210826, 0.0210332, 0.0210576, 0.0211502, 0.0212008, 0.0209904, 0.0210191, 0.0210775, 0.0211045, 0.0209879, 0.0211051, 0.0210568, 0.0209078, 0.0206866, 0.0204754, 0.0204566, 0.0204005, 0.0204629, 0.0202649, 0.0201598, 0.0202263, 0.0202866, 0.0202853, 0.0201766, 0.0201573, 0.0201199, 0.0201266, 0.0202045, 0.0203559, 0.0205059, 0.0205846, 0.020627, 0.0205889, 0.0207248, 0.0210338, 0.020848, 0.0208437, 0.0208103, 0.0208386, 0.0210152, 0.0210037, 0.0209833, 0.0209339, 0.0211464, 0.0211956, 0.021223, 0.0212042, 0.0212204, 0.0212783, 0.0211757, 0.0211111, 0.0210839, 0.0210912, 0.0210711, 0.0210395, 0.0209084, 0.0208914, 0.0208499, 0.0208768, 0.0208232, 0.0207477, 0.0206863, 0.0205907, 0.0204855, 0.0203648, 0.0202348, 0.0201752, 0.0201825, 0.0202359, 0.0197321, 0.0197388, 0.0197846, 0.019871, 0.0199399, 0.0199271, 0.0199092, 0.0199943, 0.0199323, 0.019925, 0.0199468, 0.0198981, 0.0197694, 0.0196401, 0.0196609, 0.0197674, 0.0199079, 0.0198764, 0.0197805, 0.0196566, 0.019635, 0.0196172, 0.0196139, 0.0195516, 0.0194965, 0.0194379, 0.0193647, 0.0194341, 0.0195323, 0.0194786, 0.0194654, 0.0194412, 0.0195673, 0.0196621, 0.0200041, 0.0201452, 0.0203318, 0.0207462, 0.020807, 0.0208663, 0.0211534, 0.0214383, 0.0215592, 0.0217233, 0.0218282, 0.0218163, 0.0218591, 0.0218237, 0.021759, 0.0217767, 0.0217857, 0.0218089, 0.0217847, 0.0217291, 0.021772, 0.0217833, 0.0217315, 0.0216713, 0.0215465, 0.021455, 0.0214184, 0.0213807, 0.0213498, 0.0214956, 0.0213265, 0.0211701, 0.0211106, 0.0212494, 0.0211715, 0.0211351, 0.0211105, 0.0210394, 0.0210233, 0.0210684, 0.0210343, 0.0210338, 0.0211122, 0.0211723, 0.0212455, 0.0213314, 0.021379, 0.0213257, 0.0212567, 0.021173, 0.0210044, 0.020943, 0.0209889, 0.021074, 0.0211136, 0.0211298, 0.0210848, 0.0208255, 0.0207487, 0.0205689, 0.0203249, 0.0201748, 0.0200892, 0.0199557, 0.0199712, 0.0200732, 0.0200941, 0.0200531, 0.0202368, 0.0202247, 0.0202549, 0.0204321, 0.0207082, 0.0209148, 0.0211646, 0.0214232, 0.0217762, 0.0220687, 0.0222249, 0.0224237, 0.0223453, 0.0222845, 0.0222731, 0.0223779, 0.0225421, 0.0226999, 0.0227258, 0.0226618, 0.0225859, 0.0224936, 0.0222001, 0.0221718, 0.0221996, 0.0221581, 0.0221579, 0.0222042, 0.0221501, 0.022077, 0.0219974, 0.021962, 0.0218739, 0.0217455, 0.0216513, 0.0215163, 0.0213792, 0.0213222, 0.0212691, 0.0211987, 0.0211547, 0.0210673, 0.0210157, 0.0211074, 0.0212388, 0.0213335, 0.0213757, 0.0215502, 0.0216351, 0.0217362, 0.0217988, 0.0217344, 0.0217221, 0.0216904, 0.0218072, 0.0219162, 0.0219047, 0.0218265, 0.0217598, 0.0216885, 0.0216105, 0.0214485, 0.0212917, 0.0211206, 0.0207932, 0.0207192, 0.0206767, 0.0205668, 0.0204963, 0.0204468, 0.0204232, 0.0206817, 0.0205926, 0.0206728, 0.0206639, 0.0207806, 0.0208969, 0.0210323, 0.0212097, 0.0209803, 0.0211011, 0.0212169, 0.0213352, 0.0215127, 0.0216621, 0.0217956, 0.0218957, 0.0219924, 0.0220033, 0.021993, 0.0220259, 0.0220868, 0.0222158, 0.0221891, 0.0222084, 0.0221308, 0.0219296, 0.0217949, 0.0217846, 0.0218326, 0.0219313, 0.0219141, 0.0218403, 0.0217735, 0.0211396, 0.0210943, 0.021037, 0.0209607, 0.0209538, 0.0209682, 0.020957, 0.0209123, 0.0209004, 0.0207769, 0.0206586, 0.0206507, 0.0206338, 0.0206661, 0.0207139, 0.0206755, 0.0208063, 0.0208376, 0.0209504, 0.0208902 };
const double myS_v[] = { 92.2886, 90.9589, 89.6619, 88.3985, 87.1697, 85.9764, 84.8193, 83.6991, 82.6165, 81.572, 80.5662, 79.5996, 78.6727, 77.7857, 76.9392, 76.1333, 75.3686, 74.6451, 73.9632, 73.3232, 72.7253, 72.1698, 71.6568, 71.1868, 70.7599, 70.3764, 70.0366, 69.7408, 69.4894, 69.2828, 69.1213, 69.0053, 68.9354, 68.9122, 68.936, 69.0076, 69.1276, 69.2967, 69.5157, 69.7855, 70.1069, 70.4809, 70.9086, 71.391, 71.9293, 72.5248, 73.1788, 73.8927, 74.6681, 75.5065, 76.4095, 77.3791, 78.417, 79.5251, 80.7056, 81.9606, 83.2923, 84.703, 86.1953, 87.7716, 89.4346, 91.187, 93.0317, 94.9715, 97.0096, 99.149, 101.393, 103.745, 106.207, 108.785, 111.48, 114.297, 117.239, 120.31, 123.514, 126.853, 130.333, 133.955, 137.725, 141.645, 145.718, 149.95, 154.341, 158.897, 163.618, 168.509, 173.572, 178.809, 184.222, 189.812, 195.581, 201.529, 207.657, 213.966, 220.453, 227.119, 233.961, 240.978, 248.165, 255.521, 263.039, 270.715, 278.543, 286.516, 294.626, 302.865, 311.224, 319.692, 328.259, 336.911, 345.636, 354.421, 363.25, 372.108, 380.978, 389.844, 398.688, 407.49, 416.233, 424.895, 433.458, 441.9, 450.2, 458.339, 466.293, 474.043, 481.566, 488.843, 495.853, 502.575, 508.99, 515.079, 520.824, 526.208, 531.214, 535.827, 540.033, 543.82, 547.175, 550.09, 552.555, 554.564, 556.111, 557.193, 557.807, 557.953, 557.631, 556.846, 555.6, 553.901, 551.755, 549.172, 546.162, 542.736, 538.907, 534.69, 530.1, 525.152, 519.864, 514.254, 508.34, 502.14, 495.676, 488.965, 482.028, 474.885, 467.557, 460.063, 452.424, 444.658, 436.786, 428.826, 420.796, 412.715, 404.6, 396.467, 388.332, 380.211, 372.118, 364.067, 356.071, 348.141, 340.29, 332.528, 324.864, 317.307, 309.866, 302.547, 295.359, 288.305, 281.392, 274.624, 268.006, 261.539, 255.227, 249.073, 243.077, 237.241, 231.565, 226.05, 220.696, 215.502, 210.468, 205.591, 200.872, 196.308, 191.897, 187.638, 183.527, 179.564, 175.745, 172.068, 168.53, 165.129, 161.861, 158.725, 155.716, 152.833, 150.072, 147.431, 144.906, 142.496, 140.196, 138.006, 135.921, 133.939, 132.058, 130.275, 128.589, 126.995, 125.493, 124.079, 122.752, 121.51, 120.35, 119.271, 118.271, 117.347, 116.499, 115.723, 115.02, 114.386, 113.822, 113.324, 112.891, 112.523, 112.218, 111.974, 111.791, 111.666, 111.599, 111.589, 111.634, 111.733, 111.885, 112.089, 112.344, 112.648, 113.001, 113.401, 113.847, 114.338, 114.872, 115.449, 116.067, 116.725, 117.421, 118.154, 118.923, 119.725, 120.56, 121.426, 122.321, 123.243, 124.191, 125.163, 126.156, 127.169, 128.199, 129.244, 130.303, 131.372, 132.45, 133.534, 134.62, 135.708, 136.793, 137.874, 138.948, 140.011, 141.061, 142.095, 143.11, 144.103, 145.071, 146.012, 146.922, 147.799, 148.639, 149.44, 150.199, 150.914, 151.581, 152.198, 152.763, 153.274, 153.728, 154.123, 154.458, 154.73, 154.938, 155.081, 155.156, 155.165, 155.104, 154.975, 154.775, 154.506, 154.167, 153.757, 153.279, 152.731, 152.115, 151.432, 150.682, 149.868, 148.991, 148.052, 147.054, 145.998, 144.886, 143.722, 142.507, 141.244, 139.935, 138.583, 137.192, 135.764, 134.302, 132.808, 131.287, 129.74, 128.172, 126.584, 124.981, 123.364, 121.737, 120.103, 118.465, 116.825, 115.186, 113.55, 111.921, 110.301, 108.691, 107.095, 105.514, 103.951, 102.407, 100.885, 99.3853, 97.9106, 96.4621, 95.0414, 93.6498 };

void readVectorFromFile( const char* file, vector<double>& vec ){
    ifstream str( file );
    if( !str.good() ){
        cerr << "Unable to read file " << file << endl;
        assert(false);
    }
    double val;
    while( str >> val )
        vec.push_back( val );
}

class MosqLifeCycleSuite : public CxxTest::TestSuite
{
public:
    MosqLifeCycleSuite () {
        readVectorFromFile( "MLCS_EmergenceRateConstantRA.txt", adultEmergenceConstRA );
        readVectorFromFile( "MLCS_AdultMosqConstantRA.txt", adultMosquitoesConstRA );
        readVectorFromFile( "MLCS_OvipositingAdultsConstantRA.txt", ovipositingMosquitoesConstRA );
        readVectorFromFile( "MLCS_EmergenceRatePeriodicRA.txt", adultEmergencePeriodicRA );
        readVectorFromFile( "MLCS_AdultMosqPeriodicRA.txt", adultMosquitoesPeriodicRA );
        readVectorFromFile( "MLCS_OvipositingAdultsPeriodicRA.txt", ovipositingMosquitoesPeriodicRA );
        readVectorFromFile( "MLCS_PeriodicResourceAvailability.txt", periodicResourceAvailability );
        
        assert( adultEmergenceConstRA.size() == ovipositingMosquitoesConstRA.size() );
        assert( adultEmergenceConstRA.size() == adultMosquitoesConstRA.size() );
        assert( adultEmergenceConstRA.size() == 1000 );	// can easily be changed
        assert( adultEmergencePeriodicRA.size() == ovipositingMosquitoesPeriodicRA.size() );
        assert( adultEmergencePeriodicRA.size() == adultMosquitoesPeriodicRA.size() );
        assert( adultEmergencePeriodicRA.size() == 1000 ); // can easily be changed
        assert( periodicResourceAvailability.size() == YEAR_LEN );
        
        // We only actually need 3 values from mosqElt to initialise MosquitoTransmission class
        scnXml::InputValue nanIV( numeric_limits<double>::quiet_NaN() );
        scnXml::BetaMeanSample nanBMS( numeric_limits<double>::quiet_NaN(),
                                       numeric_limits<double>::quiet_NaN() );
        scnXml::Mosq mosqElt( 3, 11, nanIV, nanIV, nanIV, nanIV,
                              nanBMS, nanBMS, nanBMS, nanIV, nanIV, 0.1 );
        
        scnXml::MosqStage eggStage( 2, 0.9 );
        scnXml::LarvalStage larvalStage( 8, 0.7 );
        scnXml::MosqStage pupalStage( 1, 0.95 );
        larvalStage.getDaily().push_back( scnXml::Daily( 0.3, 1.0 ) );
        larvalStage.getDaily().push_back( scnXml::Daily( 0.4, 0.95 ) );
        larvalStage.getDaily().push_back( scnXml::Daily( 0.5, 0.9 ) );
        larvalStage.getDaily().push_back( scnXml::Daily( 0.6, 0.85 ) );
        larvalStage.getDaily().push_back( scnXml::Daily( 0.7, 0.8 ) );
        larvalStage.getDaily().push_back( scnXml::Daily( 0.8, 0.75 ) );
        larvalStage.getDaily().push_back( scnXml::Daily( 0.9, 0.7 ) );
        larvalStage.getDaily().push_back( scnXml::Daily( 1.0, 0.65 ) );
        scnXml::InputValue femaleEggsLaid( 50.0 );
        scnXml::LifeCycle lcElt( eggStage, larvalStage, pupalStage, femaleEggsLaid );
        lcElt.setEstimatedLarvalResources( 1e5 );
        ::xsd::cxx::tree::optional<scnXml::LifeCycle> lcOpt(lcElt);
        UnittestUtil::MosqLifeCycle_init();
        mt.initialise(lcOpt, mosqElt);
        LCEmergence* lce = dynamic_cast<LCEmergence*>(mt.emergence.get());
        assert (lce != 0);
        lcp = &lce->lcParams;
        lc = &lce->lifeCycle;
    }
    ~MosqLifeCycleSuite () {
    }
    
    void setUp () {
    }
    void tearDown () {
    }
    
    void testScaleVector1to1() {
        gsl_vector *source = gsl_vector_alloc( 1 );
        gsl_vector_set( source, 0, MLCS_a );
        gsl_vector *target = gsl_vector_alloc( 1 );
        vector_scale_length( source, target );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 0 ), MLCS_a );
        gsl_vector_free( source );
        gsl_vector_free( target );
    }
    void testScaleVector1to2() {
        gsl_vector *source = gsl_vector_alloc( 1 );
        gsl_vector_set( source, 0, MLCS_a );
        gsl_vector *target = gsl_vector_alloc( 2 );
        vector_scale_length( source, target );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 0 ), MLCS_a );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 1 ), MLCS_a );
        gsl_vector_free( source );
        gsl_vector_free( target );
    }
    void testScaleVector2to1() {
        gsl_vector *source = gsl_vector_alloc( 2 );
        gsl_vector_set( source, 0, MLCS_a );
        gsl_vector_set( source, 1, MLCS_b );
        gsl_vector *target = gsl_vector_alloc( 1 );
        vector_scale_length( source, target );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 0 ), (MLCS_a+MLCS_b) / 2.0 );
        gsl_vector_free( source );
        gsl_vector_free( target );
    }
    void testScaleVector3to4() {
        gsl_vector *source = gsl_vector_alloc( 3 );
        gsl_vector_set( source, 0, MLCS_a );
        gsl_vector_set( source, 1, MLCS_b );
        gsl_vector_set( source, 2, MLCS_c );
        gsl_vector *target = gsl_vector_alloc( 4 );
        vector_scale_length( source, target );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 0 ), MLCS_a );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 1 ), (0.25*MLCS_a + 0.5*MLCS_b) / 0.75 );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 2 ), (0.5*MLCS_b + 0.25*MLCS_c) / 0.75 );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 3 ), MLCS_c );
        gsl_vector_free( source );
        gsl_vector_free( target );
    }
    void testScaleVector5to3() {
        gsl_vector *source = gsl_vector_alloc( 5 );
        gsl_vector_set( source, 0, MLCS_a );
        gsl_vector_set( source, 1, MLCS_b );
        gsl_vector_set( source, 2, MLCS_c );
        gsl_vector_set( source, 3, MLCS_d );
        gsl_vector_set( source, 4, MLCS_e );
        gsl_vector *target = gsl_vector_alloc( 2 );
        vector_scale_length( source, target );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 0 ), (MLCS_a+MLCS_b+0.5*MLCS_c) / 2.5 );
        TS_ASSERT_EQUALS( gsl_vector_get( target, 1 ), (0.5*MLCS_c+MLCS_d+MLCS_e) / 2.5 );
        gsl_vector_free( source );
        gsl_vector_free( target );
    }
    
    // These two tests only look at the egg, larvae and pupae code
    void testELPStagesConstRA (){
        // constant resource availability
        lcp->invLarvalResources.assign( 1, 1e-8 );
        
        LifeCycle lcModel;
        lcModel.init( *lcp );
        // we start with 100,000 eggs and then run the model
        lcModel.newEggs[0] = 100000.0;
        
        for( size_t d=1; d<adultEmergenceConstRA.size(); ++d ){
            double emergence = lcModel.updateEmergence( *lcp, ovipositingMosquitoesConstRA[d], d, 0 );
            TS_ASSERT_APPROX( emergence, adultEmergenceConstRA[d] );
        }
    }
    void testELPStagesPeriodicRA (){
        // annually-periodic resource availability
        lcp->invLarvalResources = periodicResourceAvailability;
        
        LifeCycle lcModel;
        lcModel.init( *lcp );
        // we start with 100,000 eggs and then run the model
        lcModel.newEggs[0] = 100000.0;
        
        for( size_t d=1; d<adultEmergencePeriodicRA.size(); ++d ){
            double emergence = lcModel.updateEmergence( *lcp, ovipositingMosquitoesPeriodicRA[d], d, (d-1)%YEAR_LEN );
            TS_ASSERT_APPROX( emergence, adultEmergencePeriodicRA[d] );
        }
    }
#if 0
    // These two tests are similar to the above, but include N_v code
    void te stLifeCycleConstRA (){
        // constant resource availability; in this case we need a vector 365 long
        lcp->invLarvalResources.assign( 365, 1e-8 );
        
        vector<double> zeros( mt.N_v_length, 0.0 );
        mt.initState( P_A, P_df, 0.0, 0.0, zeros );
        // we start with 100,000 eggs and then run the model
        lc->newEggs[0] = 100000.0;
        
        for( size_t d=1; d<adultMosquitoesConstRA.size(); ++d ){
            mt.resetTSStats();
            mt.update( d, P_A, P_df, 0.0, false, false );
            double N_v = mt.N_v[ d % mt.N_v_length ];
            TS_ASSERT_APPROX( N_v, adultMosquitoesConstRA[d] );
            TS_ASSERT_APPROX( mt.getLastN_v0(), adultEmergenceConstRA[d] );
        }
    }
    void te stLifeCyclePeriodicRA (){
        // annually-periodic resource availability
        mt.lcParams.invLarvalResources = periodicResourceAvailability;
        
        vector<double> zeros( mt.N_v_length, 0.0 );
        mt.initState( P_A, P_df, 0.0, 0.0, zeros );
        // we start with 100,000 eggs and then run the model
        mt.lifeCycle.newEggs[0] = 100000.0;
        
        for( size_t d=1; d<adultMosquitoesPeriodicRA.size(); ++d ){
            mt.resetTSStats();
            mt.update( d, P_A, P_df, 0.0, false );
            double N_v = mt.N_v[ d % mt.N_v_length ];
            TS_ASSERT_APPROX( N_v, adultMosquitoesPeriodicRA[d] );
            TS_ASSERT_APPROX( mt.getLastN_v0(), adultEmergencePeriodicRA[d] );
        }
    }
    
    void t estSimulateWithInfections (){
        /* NOTE: this is only for debugging
        // constant resource availability; in this case we need a vector 365 long
        mt.lcParams.invLarvalResources.assign( 365, 1e-8 );
        
        mt.initState( P_A, P_df, initNvFromSv, initOvFromSv, vector<double>(&myS_v[0], &myS_v[365]) );
        
        for( size_t d=1; d<400; ++d ){
            mt.update( d, P_A, P_df, myP_dif[d % 365], true );
        }
        */
    }
    
    void te stFitting (){
        try{
            vector<double> fixedP_difVec( myP_dif, myP_dif + 365 );
            vector<double> fixedS_vVec( myS_v, myS_v + 365 );
            ResourceFitter clm( mt, P_A, P_df, initNvFromSv, initOvFromSv, false );
            clm.targetS_vWithP_dif( fixedS_vVec, fixedP_difVec );
            gsl_vector *x = gsl_vector_alloc( 1 );
            gsl_vector_set_all( x, 1e8 );
            clm.sampler( x );
            clm.fit();
            gsl_vector_free( x );
            cout<<"done fitting"<<endl;
        } catch (const OM::util::traced_exception& e) {
            cerr << "Exception: " << e.what() << endl;
            cerr << e << flush;
            TS_ASSERT(false);
        } catch (const OM::util::xml_scenario_error& e) {
            cerr << "Error in scenario XML file: " << e.what() << endl;
            TS_ASSERT(false);
        }
    }
#endif
    
private:
    MosqTransmission mt;
    LifeCycleParams *lcp;
    LifeCycle *lc;
    vector<double> ovipositingMosquitoesConstRA;
    vector<double> adultEmergenceConstRA;
    vector<double> adultMosquitoesConstRA;
    vector<double> ovipositingMosquitoesPeriodicRA;
    vector<double> adultEmergencePeriodicRA;
    vector<double> adultMosquitoesPeriodicRA;
    // we use the first value to calculate state at d=1
    vector<double> periodicResourceAvailability;
};

#endif
