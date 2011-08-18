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
#include <gsl/gsl_vector.h>


namespace OM { namespace Transmission { namespace Vector {
    // this is a private function not declared in the header — but we can still test it
    void vector_scale_length( const gsl_vector *source, gsl_vector *target );
} } }
using OM::Transmission::Vector::vector_scale_length;

// just some different constants
const double MLCS_a=5.23e-5, MLCS_b=9.4e12, MLCS_c=9.32, MLCS_d=9.34243e-2, MLCS_e=141.23;

class MosqLifeCycleSuite : public CxxTest::TestSuite
{
public:
    MosqLifeCycleSuite () {
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
    
private:
};

#endif
