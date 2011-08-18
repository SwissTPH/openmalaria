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

#include "Transmission/Vector/MosquitoLifeCycle.h"
#include "schema/entomology.h"
#include "util/vectors.h"
#include "util/errors.h"
#include "util/CommandLine.h"
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


void MosqLifeCycleParams::initMosqLifeCycle( const scnXml::LifeCycle& lifeCycle ){
    // Simple constants stored in XML:
    eggStageDuration = lifeCycle.getEggStage().getDuration();
    larvalStageDuration = lifeCycle.getLarvalStage().getDuration();
    pupalStageDuration = lifeCycle.getPupalStage().getDuration();
    // we're only interested in female eggs, hence divide by 2:
    fEggsLaidByOviposit = lifeCycle.getEggsLaidByOviposit().getValue() / 2.0;
    //NOTE: store daily or whole-stage probability of survival?
    pSurvEggStage = lifeCycle.getEggStage().getSurvival();
    pSurvDayAsLarvae = std::pow( lifeCycle.getLarvalStage().getSurvival(), 1.0 / larvalStageDuration );
    pSurvPupalStage = lifeCycle.getPupalStage().getSurvival();
    
    // constants varying by larval age; probably stored directly in XML:
    larvaeResourceUsage.reserve( larvalStageDuration );
    effectCompetitionOnLarvae.reserve( larvalStageDuration );
    const scnXml::LarvalStage::DailySequence& larvDev = lifeCycle.getLarvalStage().getDaily();
    for( scnXml::LarvalStage::DailyConstIterator it = larvDev.begin(); it!=larvDev.end(); ++it ){
        larvaeResourceUsage.push_back( it->getResourceUsage() );
        effectCompetitionOnLarvae.push_back( it->getEffectCompetition() );
    }
    
    // complex derivation: annual resource availability to larvae
    //NOTE: this is set by fitLarvalResourcesFromEmergence()
    larvalResources.resize(365);
    
    debugOutput = CommandLine::option( CommandLine::DEBUG_VECTOR_FITTING );
}

// forward declarations
int CaptiveLCModel_rootfind_sampler( const gsl_vector *x, void *params, gsl_vector *f );
double CaptiveLCModel_minimise_sampler( const gsl_vector *x, void *params );

enum CaptiveLCModelFitMethod {
    minimise, find_root
};

/// Container to run life-cycle model with fixed human inputs and run fitting
/// algorithms.
struct CaptiveLCModel {
    /** Store fixed parameters.
     * 
     * @param lcP Reference to MosqLifeCycleParams. This is updated with
     * different resource availabilities as they are tested. Other parameters
     * are not changed.
     * @param Pdf Average P_df value (assumed constant)
     * @param PA Average P_A value (assumed constant)
     * @param NvL N_v_length (as in SpeciesModel class)
     * @param mRD The duration of a feeding cycle (τ)
     */
    CaptiveLCModel(
        MosqLifeCycleParams& lcP,
        double Pdf, double PA,
        size_t NvL, size_t mRD
    ) :
        lcParams( lcP ),
        P_df( Pdf ), P_A( PA ),
        N_v_length( NvL ),
        mosqRestDuration( mRD ),
        fitTarget( FT_NONE ),
        target( 0 )
    {
        //WARNING: a bad initial value here can prevent the algorithm from working!
        //TODO: maybe initial guess should go in XML?
        // Set initial_guess.
        initial_guess = gsl_vector_alloc( 1 );
        gsl_vector_set_all( initial_guess, 1e-5 );
        buf = gsl_vector_alloc( lcParams.larvalResources.size() );
    }
    ~CaptiveLCModel() {
        gsl_vector_free( initial_guess );
        gsl_vector_free( buf );
    }
    
    /** Set emergence-rate as target. */
    void targetEmergenceRate( const vector<double>& emergeRate ){
        assert( false && "need to check inputs to captive model" );
        fitTarget = FT_EMERGENCE;
        target = &emergeRate;
    }
    /** Set S_v as target. */
    void targetS_v( const vector<double>& S_v ){
        fitTarget = FT_S_V;
        target = &S_v;
    }
    
    /** Run fitting algorithms (root-finding or minimisation). */
    void fit() {
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
    
private:
    /** Copy the best result to larval resources. fit(365) must have been
     * called previously to get the correct length. */
    inline void copyBestLarvalResources (){
        copyToLarvalResources( initial_guess );
    }
    
    /** Run root-finding/minimisation algorithm.
     *
     * @param order The number of dimensions to try fitting at once. Probably
     * only 365 is valid for root-finding (since otherwise a root may not exist);
     * in any case the last time this is called must be with order 365.
     * @param method Either minimise or find_root.
     */
    void fit( size_t order, CaptiveLCModelFitMethod method ){
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
                &CaptiveLCModel_minimise_sampler,
                static_cast<void*>(this),
                initial_guess,
                step_size );
            gsl_vector_free( step_size );
        }else{
            solver = new MultidimRootFinder(
                gsl_multiroot_fsolver_hybrids,
                order,
                &CaptiveLCModel_rootfind_sampler,
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
    
    /** Run captive model. Uses a warmup length of lcParams.getTotalDuration(),
     * then a sampling duration of 365 intervals. */
    void simulate1Year()
    {
        lcModel.init( lcParams );
        size_t lastDDifferent;
        
        //TODO: are fixed P_df and P_A values adequate? I *think* so but need to check.
        //TODO: we need P_dif to get S_v which is definitely varies annually.
        // We use average P_df and P_A values ad assume these are constant (in the
        // full model they only vary dependent on human age, human heterogeneity
        // and interventions (we assume there are no interventions).
        // TODO: I think we can just initialise this from 0 and same for O_v/S_v?
        vector<double> N_v( N_v_length );
        
        size_t start_day = 0;	// historical utility but possibly still useful
        size_t sample_start = start_day + lcParams.getTotalDuration();
        size_t end = sample_start + 10*TimeStep::DAYS_IN_YEAR;
        for( size_t d = start_day; d<end; ++d ){
            size_t dMod = d + N_v_length;
            assert (dMod >= (size_t)N_v_length);
            // Indecies for today, yesterday and mosqRestDuration days back:
            size_t t    = dMod % N_v_length;
            size_t t1   = (dMod - 1) % N_v_length;
            size_t ttau = (dMod - mosqRestDuration) % N_v_length;
            // Day of year. Note that emergence during day 1
            // comes from mosqEmergeRate[0], hence subtraction by 1.
            size_t dYear1 = (d - 1) % TimeStep::DAYS_IN_YEAR;
            
            // update life-cycle model
            //NOTE: dependence on N_v ttau days ago — without proper initialisation this requires warmup.
            double newAdults = lcModel.updateEmergence( lcParams,
                                                        P_df * N_v[ttau],
                                                        d, dYear1
                                                        );
            
            // num seeking mosquitos is: new adults + those which didn't find a host
            // yesterday + those who found a host tau days ago and survived cycle:
            N_v[t] = newAdults
                    + P_A  * N_v[t1]
                    + P_df * N_v[ttau];
            
            if( fitTarget == FT_EMERGENCE ){
                samples[ dYear1 ] = newAdults;
            }else if( fitTarget == FT_S_V ){
                double s;       //FIXME: need to calculate S_v
                if( samples[ dYear1 ] - s > 0.001 * max(samples[dYear1],s) )
                    lastDDifferent = d;
                samples[ dYear1 ] = s;
                cout<<d<<'\t'<<s<<endl;
                if( d - lastDDifferent >= TimeStep::DAYS_IN_YEAR ){
                    cerr<<"completed first sample!"<<endl;
                    exit(1);//FIXME: just want to stop here
                    return;	// we're done
                }
            }else{
                assert( false );
            }
        }
        
        // if we get to here, for some reason the dynamic system never converged to a stable periodic orbit
        throw TRACED_EXCEPTION( "larvae resource fitting: system doesn't converge to a stable orbit", Error::VectorFitting );
    }
    
    /** Sample: given descriptor for resource availability x, calculate
     * resultant emergence rate.
     * 
     * Function sets buf to (sampled emergence rate) - (target emergence rate).
     * 
     * @param x Periodic descriptor of annual resource availability. If not of
     * length 365 days, vector will be scaled.
     */
    void sampler( const gsl_vector *x ){
        if( x->size == lcParams.larvalResources.size() ){
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
            msg << "; mean input was " << vectors::mean( lcParams.larvalResources );
            throw TRACED_EXCEPTION( msg.str(), Error::VectorFitting );
        }
    }
    
    /** Copy the exponential of input to larvalResources. */
    void copyToLarvalResources( const gsl_vector* input ){
        std::vector<double>& lr = lcParams.larvalResources;
        assert( input->size == lr.size() );
        //TODO: invert (1/x) to help algorithm?
        for( size_t i=0; i<lr.size(); ++i )
            lr[i] = gsl_vector_get( input, i );
    }
    
    MosqLifeCycleParams& lcParams;
    MosquitoLifeCycle lcModel;
    double P_df, P_A;
    size_t N_v_length, mosqRestDuration;
    enum FitTarget {
        FT_NONE,	// none set yet
        FT_EMERGENCE,
        FT_S_V
    } fitTarget;
    const std::vector<double> *target;
    // A vector which is filled with sampled values. Must be of length 1 year;
    // first value corresponds to emergence sampled for 1st day of year.
    std::vector<double> samples;
    gsl_vector *initial_guess;
    gsl_vector *buf;	// working memory, of length 365
    
    friend int CaptiveLCModel_rootfind_sampler( const gsl_vector *x, void *params, gsl_vector *f );
    friend double CaptiveLCModel_minimise_sampler( const gsl_vector *x, void *params );
};

// helper functions: GSL algorithms can only call global functions
int CaptiveLCModel_rootfind_sampler( const gsl_vector *x, void *params, gsl_vector *f ){
    CaptiveLCModel *clmp = static_cast<CaptiveLCModel *>( params );
    clmp->sampler( x );
    vector_scale_length( clmp->buf, f );
    return GSL_SUCCESS;	// we use exceptions to report errors, not error codes
}
double CaptiveLCModel_minimise_sampler( const gsl_vector *x, void *params ){
    CaptiveLCModel *clmp = static_cast<CaptiveLCModel *>( params );
    clmp->sampler( x );
    double sumSquares = 0.0;
    for( size_t i=0; i < clmp->buf->size; ++i ){
        double diff = gsl_vector_get( clmp->buf, i );
        sumSquares += diff * diff;
    }
    return sumSquares;
}

void MosqLifeCycleParams::fitLarvalResourcesFromS_v(
    const MosquitoLifeCycle& lcModel,
    double P_df, double P_A,
    size_t N_v_length, size_t mosqRestDuration,
    vector<double>& annualP_dif,
    vector<double>& targetS_v )
{
    CaptiveLCModel clm( *this, P_df, P_A, N_v_length, mosqRestDuration );
    clm.targetS_v( targetS_v );
    clm.fit();
}


void MosquitoLifeCycle::init( const MosqLifeCycleParams& lcParams ){
    // Shouldn't matter that values start at 0, since the outputs of this model
    // aren't used before all zeros have been replaced.
    newEggs.assign( lcParams.eggStageDuration, 0.0 );
    numLarvae.assign( lcParams.larvalStageDuration, 0.0 );
    newPupae.assign( lcParams.pupalStageDuration, 0.0 );
}

double MosquitoLifeCycle::getResRequirements( const MosqLifeCycleParams& lcParams ) const{
    double resReq = 0.0;
    for( int age=0; age<lcParams.larvalStageDuration; ++age){
        resReq += lcParams.larvaeResourceUsage[age] * numLarvae[age];
    }
    return resReq;
}

double MosquitoLifeCycle::updateEmergence( const MosqLifeCycleParams& lcParams,
                                           double nOvipositingMosqs,
                                           size_t d, size_t dYear1 ){
    // num newly emerging adults comes from num new pupae
    // pupalStageDuration days ago:
    double newAdults =
        lcParams.pSurvPupalStage * newPupae[d % lcParams.pupalStageDuration];
    
    // resource competition during last time-step (L(t) * gamma(t))
    double resourceCompetition = getResRequirements( lcParams )
        * lcParams.larvalResources[dYear1];
    // num new pupae uses larval development formula based on num larvae
    // which were one day away from becoming adults yesterday
    newPupae[d % lcParams.pupalStageDuration] =
        lcParams.pSurvDayAsLarvae * numLarvae[lcParams.larvalStageDuration-1] /
        ( 1.0 + resourceCompetition * lcParams.effectCompetitionOnLarvae[lcParams.larvalStageDuration-1] );
    for( size_t age=lcParams.larvalStageDuration-1;age>=1;--age ){
        numLarvae[age] = lcParams.pSurvDayAsLarvae * numLarvae[age-1] /
        ( 1.0 + resourceCompetition * lcParams.effectCompetitionOnLarvae[age-1] );
    }
    
    // num new larvae comes from num eggs laid eggStageDuration days ago:
    numLarvae[ 0 ] =
        lcParams.pSurvEggStage * newEggs[d % lcParams.eggStageDuration];
    
    // num eggs laid depends on number of mosquitoes which completed a
    // feeding & egg-laying cycle starting tau days ago:
    newEggs[d % lcParams.eggStageDuration] =
        lcParams.fEggsLaidByOviposit * nOvipositingMosqs;
    
    return newAdults;
}

}
}
}
