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
#include "Transmission/Vector/MosquitoLifeCycle.h"
#include "schema/entomology.h"
#include <gsl/gsl_vector.h>
#include <fstream>


namespace OM { namespace Transmission { namespace Vector {
    // this is a private function of the ResourceFitter module not declared in
    // the header — but we can still test it
    void vector_scale_length( const gsl_vector *source, gsl_vector *target );
} } }
using namespace OM::Transmission::Vector;

// just some different constants
const double MLCS_a=5.23e-5, MLCS_b=9.4e12, MLCS_c=9.32, MLCS_d=9.34243e-2, MLCS_e=141.23;

class MosqLifeCycleSuite : public CxxTest::TestSuite
{
public:
    MosqLifeCycleSuite () {
        double val;
        ifstream str( "MLCS_EmergenceRateConstantRA.txt" );
        assert( str.good() );
        while( str >> val )
            adultEmergenceConstRA.push_back( val );
        str.close(); str.clear();
        str.open( "MLCS_OvipositingAdultsConstantRA.txt" );
        assert( str.good() );
        while( str >> val )
            ovipositingMosquitoesConstRA.push_back( val );
        str.close(); str.clear();
        str.open( "MLCS_EmergenceRatePeriodicRA.txt" );
        assert( str.good() );
        while( str >> val )
            adultEmergencePeriodicRA.push_back( val );
        str.close(); str.clear();
        str.open( "MLCS_OvipositingAdultsPeriodicRA.txt" );
        assert( str.good() );
        while( str >> val )
            ovipositingMosquitoesPeriodicRA.push_back( val );
        str.close(); str.clear();
        str.open( "MLCS_PeriodicResourceAvailability.txt" );
        assert( str.good() );
        while( str >> val )
            periodicResourceAvailability.push_back( val );
        
        assert( adultEmergenceConstRA.size() == ovipositingMosquitoesConstRA.size() );
        assert( adultEmergenceConstRA.size() == 1000 );	// can easily be changed
        assert( adultEmergencePeriodicRA.size() == ovipositingMosquitoesPeriodicRA.size() );
        assert( adultEmergencePeriodicRA.size() == 1000 ); // can easily be changed
        assert( periodicResourceAvailability.size() == YEAR_LEN );
        
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
        lcParams.initMosqLifeCycle( lcElt );
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
    
    void testLifeCycleConstRA (){
        lcParams.invLarvalResources.assign( 1, 1e-8 );
        
        MosquitoLifeCycle lcModel;
        lcModel.init( lcParams );
        // we start with 100,000 eggs and then run the model
        lcModel.newEggs[0] = 100000.0;
        
        for( size_t d=1; d<adultEmergenceConstRA.size(); ++d ){
            double emergence = lcModel.updateEmergence( lcParams, ovipositingMosquitoesConstRA[d], d, 0 );
            TS_ASSERT_APPROX( emergence, adultEmergenceConstRA[d] );
        }
    }
    
    void testLifeCyclePeriodicRA (){
        lcParams.invLarvalResources = periodicResourceAvailability;
        
        MosquitoLifeCycle lcModel;
        lcModel.init( lcParams );
        // we start with 100,000 eggs and then run the model
        lcModel.newEggs[0] = 100000.0;
        
        for( size_t d=1; d<adultEmergencePeriodicRA.size(); ++d ){
            double emergence = lcModel.updateEmergence( lcParams, ovipositingMosquitoesPeriodicRA[d], d, (d-1)%YEAR_LEN );
            TS_ASSERT_APPROX( emergence, adultEmergencePeriodicRA[d] );
            cout << emergence << endl;
        }
    }
    
private:
    static const size_t YEAR_LEN = 365;
    MosqLifeCycleParams lcParams;
    vector<double> ovipositingMosquitoesConstRA;
    vector<double> adultEmergenceConstRA;
    vector<double> ovipositingMosquitoesPeriodicRA;
    vector<double> adultEmergencePeriodicRA;
    // we use the first value to calculate state at d=1
    vector<double> periodicResourceAvailability;
};

#endif
