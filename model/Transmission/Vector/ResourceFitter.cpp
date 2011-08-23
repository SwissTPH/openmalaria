/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "Transmission/Vector/ResourceFitter.h"
#include "util/CommandLine.h"
#include "util/vectors.h"
#include "util/errors.h"
#include "util/MultidimSolver.h"
#include <cmath>
#include <sstream>

namespace OM {
namespace Transmission {
namespace Vector {
using namespace OM::util;

bool debugOutput = false;

/** Assume slots in arrays correspond to interval [ i/l, (i+1)/l ) for slot i,
 * where l is the length of the vector. Calculates the range each slot in the 
 * target array corresponds to in the source and integrates over this range.
 * 
 * Sum of target will not equal sum of source, but the mean value should be the
 * same. **/
void vector_scale_length( const gsl_vector *source, gsl_vector *target ){
    double sf = source->size / static_cast<double>( target->size );
    for( size_t ti=0; ti<target->size; ++ti ){
        // Exact end-points in source of target cell:
        double start = sf * ti;
        double end = sf * (ti + 1);
        // Indecies in source corresponding to start and end of target cell:
        // (largest possible iEnd is source->size)
        size_t iStart = std::floor(start), iEnd = std::floor(end);
        if( iStart == iEnd ){
            // Target cell corresponds to one source cell: take that value.
            double val = gsl_vector_get( source, iStart );
            gsl_vector_set( target, ti, val );
        }else /* iStart < iEnd */{
            // Target cell corresponds to two or more souce cells: take the
            // weighted sum of all source cells.
            // Calculate weights:
            double wStart = std::floor(start) + 1.0 - start;
            double wEnd = end - std::floor(end);
            double wSum = wStart + wEnd + (iEnd - iStart - 1);
            // Note: iEnd *may* be source->size; all other indecies will be smaller
            double wVal = wStart*gsl_vector_get(source, iStart) +
                wEnd*gsl_vector_get(source, iEnd % source->size);
            for( size_t i = iStart+1; i<iEnd; ++i ){
                wVal += gsl_vector_get(source, i);
            }
            gsl_vector_set( target, ti, wVal / wSum );
        }
    }
}

// helper functions: GSL algorithms can only call global functions
int ResourceFitter_rootfind_sampler( const gsl_vector *x, void *params, gsl_vector *f ){
    ResourceFitter *clmp = static_cast<ResourceFitter *>( params );
    clmp->sampler( x );
    vector_scale_length( clmp->buf, f );
    return GSL_SUCCESS; // we use exceptions to report errors, not error codes
}
double ResourceFitter_minimise_sampler( const gsl_vector *x, void *params ){
    ResourceFitter *clmp = static_cast<ResourceFitter *>( params );
    clmp->sampler( x );
    double sumSquares = 0.0;
    for( size_t i=0; i < clmp->buf->size; ++i ){
        double diff = gsl_vector_get( clmp->buf, i );
        sumSquares += diff * diff;
    }
    return sumSquares;
}

ResourceFitter::ResourceFitter( MosquitoTransmission mosqTrans,
                                double PA, double Pdf,
                                double iNvSv, double iOvSv
                              ) :
    larvalResources( mosqTrans.getLCParams().larvalResources ),
    mosquitoTransmission( mosqTrans ),
    P_A( PA ), P_df( Pdf ),
    initNvFromSv( iNvSv ), initOvFromSv( iOvSv ),
    fitTarget( FT_NONE ),
    target( 0 )
{
    //WARNING: a bad initial value here can prevent the algorithm from working!
    //TODO: maybe initial guess should go in XML?
    // Set initial_guess.
    initial_guess = gsl_vector_alloc( 1 );
    gsl_vector_set_all( initial_guess, 1e-5 );
    buf = gsl_vector_alloc( larvalResources.size() );
    
    debugOutput = CommandLine::option( CommandLine::DEBUG_VECTOR_FITTING );
}
ResourceFitter::~ResourceFitter() {
    gsl_vector_free( initial_guess );
    gsl_vector_free( buf );
}

double assertSimilarP_dif( double avg, double x, double tol ){
    double xa = x/avg;
    if( !( (1.0/tol) < xa && xa < tol ) ){
        cerr<<"avg: "<<avg<<"; x: "<<x<<"; tol: "<<tol<<endl;
        throw TRACED_EXCEPTION( "P_dif has not converged to a fixed annual periodic cycle during initialisation", Error::VectorWarmup );
    }
    return max( xa/tol, tol/xa );
}
void ResourceFitter::targetS_vWithP_dif( const vector<double>& S_v,
                                         const vector<double>& sampledP_dif ){
    fitTarget = FT_S_V;
    target = &S_v;
    
    annualP_dif.resize( TimeStep::intervalsPerYear.asInt() );
    assert( sampledP_dif.size() % annualP_dif.size() == 0 );
    for( size_t i=0; i<sampledP_dif.size(); ++i )
        annualP_dif[ i%annualP_dif.size() ] += sampledP_dif[ i ];
    double factor = static_cast<double>(annualP_dif.size()) / static_cast<double>(sampledP_dif.size());
    for( size_t i=0; i<annualP_dif.size(); ++i )
        annualP_dif[ i ] *= factor;
    
    vector<double> sumAnnualP_dif( sampledP_dif.size() / annualP_dif.size() );
    double maxTolNeeded = 0.0;
    for( size_t i=0; i<sampledP_dif.size(); ++i ){
        //TODO: increase tolerance
        double thisTol = assertSimilarP_dif( annualP_dif[ i%annualP_dif.size() ], sampledP_dif[ i ], 2 );
        maxTolNeeded = max( maxTolNeeded, thisTol );
        sumAnnualP_dif[ i/annualP_dif.size() ] += sampledP_dif[ i ];
    }
    for( size_t y=1; y< sumAnnualP_dif.size(); ++y )
        assertSimilarP_dif( sumAnnualP_dif[0], sumAnnualP_dif[y], 2 );
}


void ResourceFitter::targetEmergenceRate( const vector<double>& emergeRate ){
    assert( false && "need to check inputs to captive model" );
    fitTarget = FT_EMERGENCE;
    target = &emergeRate;
}

void ResourceFitter::fit() {
    //TODO: switch to fitting logarithm of values so we can guarantee we don't get negative resources
    // But this doesn't appear to fit!
    //TODO: find out which exceptions this throws when it can't fit, catch them, and use minimisation approach instead.
    {
        fit( 365, find_root );
    }
    fit( 1, minimise );
    fit( 3, minimise );
    fit( 10, minimise );
    fit( 34, minimise );
    fit( 112, minimise );
    fit( 365, minimise );
    copyBestLarvalResources();
}

void ResourceFitter::fit( size_t order, FitMethod method ){
    assert( method == minimise || method == find_root );
    
    if( order != initial_guess->size ){
        gsl_vector *old_estimate = initial_guess;
        initial_guess = gsl_vector_alloc ( order );
        vector_scale_length( old_estimate, initial_guess );
        gsl_vector_free( old_estimate );
    }
    
    /* Initialize method and iterate */
    MultidimSolver *solver;
    if( method == minimise ){
        //NOTE: I don't know how best to set this parameter. It's not well documented.
        gsl_vector *step_size = gsl_vector_alloc( order );
        gsl_vector_set_all( step_size, 1.0 );
        solver = new MultidimMinimiser(
            gsl_multimin_fminimizer_nmsimplex2,
            order,
            &ResourceFitter_minimise_sampler,
            static_cast<void*>(this),
            initial_guess,
            step_size );
        gsl_vector_free( step_size );
    }else{
        solver = new MultidimRootFinder(
            gsl_multiroot_fsolver_hybrids,
            order,
            &ResourceFitter_rootfind_sampler,
            static_cast<void*>(this),
            initial_guess );
    }
    
    size_t iter;
    enum FitStatus {
        in_progress, cant_improve, success
    } fit_status = in_progress;
    for( iter=0; iter<5000000; ++iter ) {
        int status = solver->iterate();
        if (status) {
            if( status==GSL_ENOPROG ){
                fit_status = cant_improve;
                break;
            }else{
                ostringstream msg;
                msg << "[while fitting vector parameter] " << gsl_strerror( status );
                throw TRACED_EXCEPTION( msg.str(), util::Error::GSL );
            }
        }
        
        if( solver->success( 1e-6 ) ){
            fit_status = success;
            break;
        }
    }
    
    // copy our best estimate to initial_guess ready for use next time fit is called
    gsl_vector *x = solver->x();
    gsl_vector_memcpy( initial_guess, x );
    
    if( fit_status != success ){
        ostringstream msg;
        msg << "Fitting with order "<<order<<" failed after "<<iter
            <<" steps. Mean value: "<<vectors::mean( x )
            <<". Reason: ";
        if( fit_status==cant_improve )
            msg << "can't improve";
        else
            msg << "too many iterations";
        throw TRACED_EXCEPTION( msg.str(), Error::VectorFitting );
    }
    if( debugOutput ){
        cerr << "Fitting with order "<<order<<" succeeded in "<<iter
            <<" steps. Mean value: "<<vectors::mean( x )
            <<"\nFull output: "<<gsl_vector_get( x, 0 );
        for( size_t i=1; i < x->size; ++i ){
            cerr << ','<<gsl_vector_get( x, i );
        }
        cerr<<endl;
        double sumSquares = 0.0;
        for( size_t i=0; i < buf->size; ++i ){
            double diff = gsl_vector_get( buf, i );
            sumSquares += diff * diff;
        }
        cerr<<"Measure of fit: "<<sumSquares<<endl;
    }
}

void ResourceFitter::simulate1Year()
{
    // reset state data so one run can't influence another
    assert( fitTarget == FT_S_V );      // because we use target as forcedS_v below:
    mosquitoTransmission.initState( P_A, P_df, initNvFromSv, initOvFromSv, *target );
    
    //TODO: we need P_dif to get S_v which is definitely varies annually.
    
    size_t end = 10*TimeStep::DAYS_IN_YEAR;
    size_t lastDDifferent;
    for( size_t d = 1; d<end; ++d ){
        double S_v = mosquitoTransmission.update( d, P_A, P_df, annualP_dif[d%TimeStep::DAYS_IN_YEAR] );
        
        size_t dYear = d % TimeStep::DAYS_IN_YEAR;
        if( fitTarget == FT_EMERGENCE ){
            assert(false);
            // NOTE: should this have been (d-1)%daysinyear or just d%...?
            //samples[ dYear1 ] = newAdults;
        }else if( fitTarget == FT_S_V ){
            if( samples[ dYear ] - S_v > 0.001 * max(samples[dYear], S_v) )
                lastDDifferent = d;
            samples[ dYear ] = S_v;
            cout<<d<<'\t'<<S_v<<endl;
            if( d - lastDDifferent >= TimeStep::DAYS_IN_YEAR ){
                return;     // we're done
            }
        }else{
            assert( false );	// fitTarget wasn't set?
        }
    }
    
    // if we get to here, for some reason the dynamic system never converged to a stable periodic orbit
    throw TRACED_EXCEPTION( "larvae resource fitting: system doesn't converge to a stable orbit", Error::VectorFitting );
}

void ResourceFitter::sampler( const gsl_vector *x ){
    for( size_t i=0; i<x->size; ++i )
        assert( finite( gsl_vector_get( x, i ) ) );
    if( x->size == larvalResources.size() ){
        copyToLarvalResources( x );
    }else{
        vector_scale_length( x, buf );
        copyToLarvalResources( buf );
    }
    
    //TODO: add option to fit shape of N_v instead of N_v0
    //TODO: do we want to add an option to fit the shape of O_v/S_v/EIR? If
    // so we need a human-infectiousness (kappa) input; perhaps this can be
    // constant or sampled from humans while they are being exposed to
    // forced EIR.
    samples.resize( TimeStep::DAYS_IN_YEAR );
    simulate1Year();
    
    double res = 0.0;
    for( size_t i=0; i < buf->size; ++i ){
        double diff = samples[i] - (*target)[i];
        gsl_vector_set( buf, i, diff );
        res += diff*diff;
    }
    if( debugOutput ){
        cerr << "Iteration has mean input "<<vectors::mean( x )<<"; fitness "<<res<<endl;
    }
    
    if( !finite(res) ){
        ostringstream msg;
        msg << "non-finite output with mean " << vectors::mean( samples );
        msg << "; mean input was " << vectors::mean( larvalResources );
        throw TRACED_EXCEPTION( msg.str(), Error::VectorFitting );
    }
}

void ResourceFitter::copyToLarvalResources( const gsl_vector* input ){
    assert( input->size == larvalResources.size() );
    //TODO: invert (1/x) to help algorithm?
    for( size_t i=0; i<larvalResources.size(); ++i )
        larvalResources[i] = gsl_vector_get( input, i );
}
    
}
}
}
