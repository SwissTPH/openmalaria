/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include <cmath>
#include <fstream>
#include "gsl_multimin.h"

namespace OM {
namespace Transmission {
using namespace OM::util;

void MosqLifeCycleParams::initMosqLifeCycle( const scnXml::LifeCycle& lifeCycle ){
    // Simple constants stored in XML:
    eggStageDuration = lifeCycle.getEggStage().getDuration();
    larvalStageDuration = lifeCycle.getLarvalStage().getDuration();
    pupalStageDuration = lifeCycle.getPupalStage().getDuration();
    // we're only interested in female eggs, hence divide by 2:
    fEggsLaidByOviposit = lifeCycle.getEggsLaidByOviposit().getValue() / 2.0;
    //NOTE: store daily or whole-stage probability of survival?
    pSurvEggStage = lifeCycle.getEggStage().getSurvival();
    pSurvDayAsLarvae = pow( lifeCycle.getLarvalStage().getSurvival(), 1.0 / larvalStageDuration );
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
}

// forward declaration
double CaptiveLCModel_sampler_wrapper( const gsl_vector *x, void *params );

// Container to run life-cycle model with fixed human inputs
struct CaptiveLCModel {
    /** Store fixed parameters.
     * 
     * @param lcP Reference to MosqLifeCycleParams. This is updated with
     * different resource availabilities as they are tested. Other parameters
     * are not changed.
     * @param Pdf Average P_df value (assumed constant)
     * @param PA Average P_A value (assumed constant)
     * @param sd The day (d parameter in VectorAnopheles::advancePeriod) at
     * which the next update should take place.
     * @param mRD The duration of a feeding cycle (τ)
     * @param N_v_orig N_v array from VectorAnopheles (will be copied)
     * @param mER Fixed emergence rate for forcing model and fitting target. Read-only.
     */
    CaptiveLCModel(
        MosqLifeCycleParams& lcP,
        double Pdf, double PA,
        size_t sd, size_t mRD,
        const std::vector<double>& N_v_orig,
        const std::vector<double>& mER
    ) :
        lcParams( lcP ),
        P_df( Pdf ), P_A( PA ),
        start_day( sd ), mosqRestDuration( mRD ),
        init_N_v( N_v_orig ), mosqEmergeRate( mER )
    {
        fout.open("MLC-fitting.csv");
        fout<<"target\t"<<mosqEmergeRate<<'\n';
    }
    
    /** Run captive model. Uses a warmup length of lcParams.getTotalDuration(),
     * then a sampling duration of 365 intervals. */
    void simulate1Year()
    {
        lcModel.init( lcParams );
        
        // We use average P_df and P_A values ad assume these are constant (in the
        // full model they only vary dependent on human age, human heterogeneity
        // and interventions (we assume there are no interventions).
        size_t N_v_length = init_N_v.size();
        // We take a copy of N_v.
        vector<double> N_v = init_N_v;
        
        size_t sample_start = start_day + lcParams.getTotalDuration();
        size_t end = sample_start + TimeStep::fromYears(1).inDays();
        for( size_t d = start_day; d<end; ++d ){
            size_t dMod = d + N_v_length;
            assert (dMod >= (size_t)N_v_length);
            // Indecies for today, yesterday and mosqRestDuration days back:
            size_t t    = dMod % N_v_length;
            size_t t1   = (dMod - 1) % N_v_length;
            size_t ttau = (dMod - mosqRestDuration) % N_v_length;
            // Day of year. Note that emergence during day 1
            // comes from mosqEmergeRate[0], hence subtraction by 1.
            size_t dYear1 = (d - 1) % TimeStep::fromYears(1).inDays();
            
            // update life-cycle model
            double newAdults = lcModel.updateEmergence( lcParams,
                                                        P_df * N_v[ttau],
                                                        d, dYear1
                                                        );
            sampledEmergence[ dYear1 ] = newAdults;
            
            // num seeking mosquitos is: new adults + those which didn't find a host
            // yesterday + those who found a host tau days ago and survived cycle:
            // Note: we fit with forced emergence to avoid a fully-dynamic model
            N_v[t] = mosqEmergeRate[dYear1]
                    + P_A  * N_v[t1]
                    + P_df * N_v[ttau];
        }
    }
    
    /** Sample: see how close the new resource rate regenerates our emergence
     * rate.
     */
    double sampler( const gsl_vector *x ){
        double scale_reduction = x->size / static_cast<double>( lcParams.larvalResources.size() );
        for( size_t i=0,len=lcParams.larvalResources.size(); i<len; ++i ){
            size_t index = static_cast<size_t>( i * scale_reduction );
            lcParams.larvalResources[i] = gsl_vector_get(x, index);
        }
        vectors::gsl2std( x, lcParams.larvalResources );
        // interpret x as a fourier series instead?
        
        //TODO: fix this function. Doesn't appear to return sensible values/output is always the same.
        //TODO: add option to fit shape of N_v instead of N_v0
        //TODO: do we want to add an option to fit the shape of O_v/S_v/EIR? If
        // so we need a human-infectiousness (kappa) input; perhaps this can be
        // constant or sampled from humans while they are being exposed to
        // forced EIR.
        sampledEmergence.resize( mosqEmergeRate.size() );
        simulate1Year();
        double sumSquares = 0.0;
        for( size_t i=0; i<sampledEmergence.size(); ++i ){
            double diff = sampledEmergence[i] - mosqEmergeRate[i];
            sumSquares += diff*diff;
        }
        return sumSquares;
    }
    
    /** Run root-finding/minimisation algorithm.
     *
     * @param order The number of dimensions to try fitting at once. Must
     * eventually be done with 365, but it may be faster to start with fewer.
     */
    void fit( int order ){
        //TODO: calculate a sensible starting point. Currently I have no idea what to use, but maybe this can be some function of mosqEmergeRate?
        gsl_vector *x = gsl_vector_alloc ( order );
        gsl_vector_set_all (x, 0.01);
        
        // Set initial step size: it seems something like 1/200th of this is added/subtracted each step
        gsl_vector *step_size = gsl_vector_alloc (order);
        gsl_vector_set_all (step_size, 1);
        
        /* Initialize method and iterate */
        gsl_multimin_function minex_func;
        minex_func.n = order;
        minex_func.f = &CaptiveLCModel_sampler_wrapper;
        minex_func.params = static_cast<void*>(this);
        
        gsl_multimin_fminimizer *s =
            gsl_multimin_fminimizer_alloc (gsl_multimin_fminimizer_nmsimplex2, order);
        gsl_multimin_fminimizer_set (s, &minex_func, x, step_size);
        
        size_t iter;
        bool success = false;
        for( iter=0; iter<1000; ++iter ) {
            int status = gsl_multimin_fminimizer_iterate(s);
            if (status) {
                if( status==GSL_ENOPROG ){
                    cout << "stopping iterations: can't improve" << endl;
                    break;
                }else{
                    ostringstream msg;
                    msg << "stopping with error code " << gsl_strerror( status );
                    throw TRACED_EXCEPTION( msg.str(), util::Error::GSL );
                }
            }
            
            double size = gsl_multimin_fminimizer_size (s);
            /*
            printf ("%5lu %10.3e %10.3e f() = %7.3f size = %.3f\n", 
                    iter,
                    gsl_vector_get (s->x, 0), 
                    gsl_vector_get (s->x, 1), 
                    s->fval, size);
            */
            if( gsl_multimin_test_size (size, 1e-5) == GSL_SUCCESS){
                success = true;
                break;
            }
        }
        
        gsl_vector_free(x);
        gsl_vector_free(step_size);
        gsl_multimin_fminimizer_free (s);
        
        cout<<"Fitting with dimension "<<order<<" in "<<iter<<" steps: "
            <<(success ? "success" : "failure")
            <<" (least squares: "<<s->fval<<")"<<endl;
        fout<<iter<<'\t'<<s->fval<<'\t'<<sampledEmergence<<endl;
    }
    
private:
    MosqLifeCycleParams& lcParams;
    MosquitoLifeCycle lcModel;
    double P_df, P_A;
    size_t start_day, mosqRestDuration;
    const std::vector<double>& init_N_v, mosqEmergeRate;
    // A vector which is filled with sampled emergence values. Must be of
    // length 1 year; first value corresponds to emergence sampled for 1st day
    // of year.
    std::vector<double> sampledEmergence;
    ofstream fout;
};

// helper function: GSL algorithms can only call global functions
double CaptiveLCModel_sampler_wrapper( const gsl_vector *x, void *params ){
    CaptiveLCModel *clmp = static_cast<CaptiveLCModel *>( params );
    return clmp->sampler( x );
}

void MosqLifeCycleParams::fitLarvalResourcesFromEmergence(
    const MosquitoLifeCycle& lcModel,
    double P_df, double P_A,
    size_t start_day, size_t mosqRestDuration,
    const std::vector<double>& N_v_orig,
    const std::vector<double>& mosqEmergeRate )
{
    CaptiveLCModel clm( *this, P_df, P_A, start_day, mosqRestDuration,
                    N_v_orig, mosqEmergeRate );
    
    //TODO: fit using gsl multi-minimiser
    //TODO: maybe it's sensible to reduce the order of the system somehow — with a lower-fidelity curve or fourier series?
    //TODO: switch to fitting logarithm of values so we can guarantee we don't get negative resources
    
    clm.fit( 1 );
    clm.fit( 5 );
    clm.fit( 73 );
    clm.fit( 365 );
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
