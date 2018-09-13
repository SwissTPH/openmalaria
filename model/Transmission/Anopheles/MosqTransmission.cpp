/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

#include "Transmission/Anopheles/MosqTransmission.h"
#include "Transmission/Anopheles/FixedEmergence.h"
#include "Transmission/Anopheles/SimpleMPDEmergence.h"
#include "WithinHost/Genotypes.h"
#include "mon/reporting.h"
#include "util/vectors.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include "schema/entomology.h"

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;
using WithinHost::Genotypes;

// -----  Initialisation of model, done before human warmup  ------

MosqTransmission::MosqTransmission() :
        mosqRestDuration(SimTime::zero()),
        EIPDuration(SimTime::zero()),
        N_v_length(SimTime::zero()),
        minInfectedThreshold( std::numeric_limits< double >::quiet_NaN() ),     // requires config
        timeStep_N_v0(0.0)
{
    // Warning: don't allocate memory here. The whole instance will be
    // bit-copied (see species.resize ... line in VectorModel constructor), so
    // each copy would point to the same object.
}


void MosqTransmission::initialise ( const scnXml::AnophelesParams::LifeCycleOptional& lcOpt,
                                    const scnXml::AnophelesParams::SimpleMPDOptional& simpleMPDOpt,
                                    const scnXml::Mosq& mosq ) {
    if (util::ModelOptions::option( util::VECTOR_LIFE_CYCLE_MODEL )){
        throw util::xml_scenario_error("VECTOR_LIFE_CYCLE_MODEL not yet "
            "implemented. Use VECTOR_SIMPLE_MPD_MODEL instead.");
        /*TODO
         * Note: this model is older than SimpleMPD and more complicated.
         * Difficulties are in parameterisation and estimation of resources.
        if (!lcOpt.present())
            throw util::xml_scenario_error(
                "VECTOR_LIFE_CYCLE_MODEL: requires <lifeCycle> element with "
                "model parameters for each anopheles species");
        emergence = unique_ptr<EmergenceModel>( new LCEmergence() );
        emergence->initLifeCycle( lcOpt.get() );
        */
    }else if (util::ModelOptions::option( util::VECTOR_SIMPLE_MPD_MODEL )){
        if (!simpleMPDOpt.present())
            throw util::xml_scenario_error(
                "VECTOR_SIMPLE_MPD_MODEL: requires <simpleMPD> element with "
                "model parameters for each anopheles species");
        emergence = unique_ptr<EmergenceModel>(new SimpleMPDEmergence(simpleMPDOpt.get()) );
    }else
        emergence = unique_ptr<EmergenceModel>( new FixedEmergence() );
    
    
    // -----  Set model variables  -----

    mosqRestDuration = SimTime::fromDays(mosq.getMosqRestDuration().getValue());
    EIPDuration = SimTime::fromDays(mosq.getExtrinsicIncubationPeriod().getValue());
    if (SimTime::oneDay() > mosqRestDuration || mosqRestDuration * 2 >= EIPDuration) {
        //TODO: limit was EIPDuration >= mosqRestDuration >= 1
        // but in usage of ftauArray this wasn't enough. Check why.
        throw util::xml_scenario_error ("Code expects EIPDuration > 2*mosqRestDuration >= 2");
    }
    N_v_length = EIPDuration + mosqRestDuration;
    
    minInfectedThreshold = mosq.getMinInfectedThreshold();
    
    
    // -----  allocate memory  -----
    // Set up fArray and ftauArray. Each step, all elements not set here are
    // calculated, even if they aren't directly used in the end;
    // however all calculated values are used in calculating the next value.
    fArray.resize(EIPDuration-mosqRestDuration+SimTime::oneDay());
    fArray[SimTime::zero()] = 1.0;
    ftauArray.resize(EIPDuration);
    for( SimTime i = SimTime::zero(); i < mosqRestDuration; i += SimTime::oneDay() ){
        ftauArray[i] = 0.0;
    }
    ftauArray[mosqRestDuration] = 1.0;
    uninfected_v.resize(N_v_length);
    uninfected_v[SimTime::zero()] = numeric_limits<double>::quiet_NaN();    // index not used
}

void MosqTransmission::initIterateScale ( double factor ){
    vectors::scale (N_v, factor);
    // What factor exactly these should be scaled by isn't obvious; in any case
    // they should reach stable values quickly.
    vectors::scale (O_v, factor);
    vectors::scale (S_v, factor);
}

void MosqTransmission::initState ( double tsP_A, double tsP_df, double tsP_dff,
                                       double initNvFromSv, double initOvFromSv,
                                       const vecDay<double>& forcedS_v ){
    N_v  .assign (N_v_length, numeric_limits<double>::quiet_NaN());
    O_v  .assign (N_v_length, Genotypes::N(), numeric_limits<double>::quiet_NaN());
    S_v  .assign (N_v_length, Genotypes::N(), numeric_limits<double>::quiet_NaN());
    P_A  .assign (N_v_length, numeric_limits<double>::quiet_NaN());
    P_df .assign (N_v_length, numeric_limits<double>::quiet_NaN());
    P_dif.assign (N_v_length, Genotypes::N(), 0.0);// humans start off with no infectiousness.. so just wait
    P_dff.assign (N_v_length, numeric_limits<double>::quiet_NaN());
    
    // Initialize per-day variables; S_v, N_v and O_v are only estimated
    assert( N_v_length <= forcedS_v.size() );
    for( SimTime t = SimTime::zero(); t < N_v_length; t += SimTime::oneDay() ){
        P_A[t] = tsP_A;
        P_df[t] = tsP_df;
        P_dff[t] = tsP_dff;
        N_v[t] = forcedS_v[t] * initNvFromSv;
        for( size_t genotype = 0; genotype < Genotypes::N(); ++genotype ){
            S_v.at(t, genotype) = forcedS_v[t] * Genotypes::initialFreq(genotype);
            O_v.at(t,genotype) = S_v.at(t,genotype) * initOvFromSv;
        }
    }
}


void MosqTransmission::update( SimTime d0, double tsP_A, double tsP_df,
        const vector<double> tsP_dif, double tsP_dff,
        bool isDynamic,
        vector<double>& partialEIR, double EIR_factor )
{
    SimTime d1 = d0 + SimTime::oneDay();    // end of step
    
    // We add N_v_length so that we can use mod_nn() instead of mod().
    SimTime d1Mod = d1 + N_v_length;
    assert (d1Mod >= N_v_length);
    // Indecies for end time, start time, and mosqRestDuration days before end time:
    SimTime t1    = mod_nn(d1, N_v_length);
    SimTime t0   = mod_nn(d0, N_v_length);
    SimTime ttau = mod_nn(d1Mod - mosqRestDuration, N_v_length);
    
    // These only need to be calculated once per time step, but should be
    // present in each of the previous N_v_length - 1 positions of arrays.
    P_A[t1] = tsP_A;
    P_df[t1] = tsP_df;
    P_dff[t1] = tsP_dff;
    for( size_t i = 0; i < Genotypes::N(); ++i )
        P_dif.at(t1,i) = tsP_dif[i];
    
    
    //BEGIN cache calculation: fArray, ftauArray, uninfected_v
    // Set up array with n in 1..θ_s−τ for f(d1Mod-n) (NDEMD eq. 1.6)
    for( SimTime n = SimTime::oneDay(); n <= mosqRestDuration; n += SimTime::oneDay() ){
        const SimTime tn = mod_nn(d1Mod-n, N_v_length);
        fArray[n] = fArray[n-SimTime::oneDay()] * P_A[tn];
    }
    fArray[mosqRestDuration] += P_df[ttau];
    
    const SimTime fAEnd = EIPDuration-mosqRestDuration;
    for( SimTime n = mosqRestDuration+SimTime::oneDay(); n <= fAEnd; n += SimTime::oneDay() ){
        const SimTime tn = mod_nn(d1Mod-n, N_v_length);
        fArray[n] =
            P_df[tn] * fArray[n - mosqRestDuration]
            + P_A[tn] * fArray[n-SimTime::oneDay()];
    }
    
    // Set up array with n in 1..θ_s−1 for f_τ(d1Mod-n) (NDEMD eq. 1.7)
    const SimTime fProdEnd = mosqRestDuration * 2;
    for( SimTime n = mosqRestDuration+SimTime::oneDay(); n <= fProdEnd; n += SimTime::oneDay() ){
        SimTime tn = mod_nn(d1Mod-n, N_v_length);
        ftauArray[n] = ftauArray[n-SimTime::oneDay()] * P_A[tn];
    }
    ftauArray[fProdEnd] += P_df[mod_nn(d1Mod-fProdEnd, N_v_length)];

    for( SimTime n = fProdEnd+SimTime::oneDay(); n < EIPDuration; n += SimTime::oneDay() ){
        SimTime tn = mod_nn(d1Mod-n, N_v_length);
        ftauArray[n] =
            P_df[tn] * ftauArray[n - mosqRestDuration]
            + P_A[tn] * ftauArray[n-SimTime::oneDay()];
    }
    
    for( SimTime d = SimTime::oneDay(); d < N_v_length; d += SimTime::oneDay() ){
        SimTime t = mod_nn(d1Mod - d, N_v_length);
        double sum = N_v[t];
        for( size_t i = 0; i < Genotypes::N(); ++i ) sum -= O_v.at(t,i);
        uninfected_v[d] = sum;
    }
    //END cache calculation: fArray, ftauArray, uninfected_v
    
    double total_S_v = 0.0;
    for( size_t genotype = 0; genotype < Genotypes::N(); ++genotype ){
        // Num infected seeking mosquitoes is the new ones (those who were
        // uninfected tau days ago, started a feeding cycle then, survived and
        // got infected) + those who didn't find a host yesterday + those who
        // found a host tau days ago and survived a feeding cycle.
        O_v.at(t1,genotype) = P_dif.at(ttau,genotype) * uninfected_v[mosqRestDuration]
                    + P_A[t0]  * O_v.at(t0,genotype)
                    + P_df[ttau] * O_v.at(ttau,genotype);
        
        //BEGIN S_v
        double sum = 0.0;
        const SimTime ts = d1Mod - EIPDuration;
        for( SimTime l = SimTime::oneDay(); l < mosqRestDuration; l += SimTime::oneDay() ){
            const SimTime tsl = mod_nn(ts - l, N_v_length); // index d1Mod - theta_s - l
            sum += P_dif.at(tsl,genotype) * P_df[ttau] * (uninfected_v[EIPDuration+l]) *
                    ftauArray[EIPDuration+l-mosqRestDuration];
        }
        
        const SimTime tsm = mod_nn(ts, N_v_length);       // index d1Mod - theta_s
        S_v.at(t1,genotype) = P_dif.at(tsm,genotype) *
                fArray[EIPDuration-mosqRestDuration] * (uninfected_v[EIPDuration])
            + sum
            + P_A[t0]*S_v.at(t0,genotype)
            + P_df[ttau]*S_v.at(ttau,genotype);


        if( isDynamic ){
            // We cut-off transmission when no more than X mosquitos are infected to
            // allow true elimination in simulations. Unfortunately, it may cause problems with
            // trying to simulate extremely low transmission, such as an R_0 case.
            if ( S_v.at(t1,genotype) <= minInfectedThreshold ) { // infectious mosquito cut-off
                S_v.at(t1,genotype) = 0.0;
                /* Note: could report; these reports often occur too frequently, however
                if( S_v[t1] != 0.0 ){        // potentially reduce reporting
            cerr << sim::ts0() <<":\t S_v cut-off"<<endl;
                } */
            }
        }
        
        partialEIR[genotype] += S_v.at(t1, genotype) * EIR_factor;
        total_S_v += S_v.at(t1, genotype);
        //END S_v
    }
    
    const double nOvipositing = P_dff[ttau] * N_v[ttau];       // number ovipositing on this step
    const double newAdults = emergence->update( d0, nOvipositing, total_S_v );
    util::streamValidate( newAdults );
    
    // num seeking mosquitos is: new adults + those which didn't find a host
    // yesterday + those who found a host tau days ago and survived cycle:
    N_v[t1] = newAdults
                + P_A[t0]  * N_v[t0]
                + nOvipositing;
    
    timeStep_N_v0 += newAdults;
    
//     if( printDebug ){
//         cerr<<"step ending "<<d1<<" (days):\temergence "<<newAdults<<",\tN_v "<<N_v[t1]<<",\tS_v "<<total_S_v<<endl;
        //cerr << "len: "<<N_v_length<<"\td1Mod: "<<d1Mod<<"\tt(0,1): "<<t0<<" "<<t1<<" "<<ttau<<endl;
/*        cerr<<"P_A\t"<<P_A[t0]<<"\t"<<P_A[t1]<<"\t"<<P_A[ttau]<<endl;
        cerr<<"P_df\t"<<P_df[t0]<<"\t"<<P_df[t1]<<"\t"<<P_df[ttau]<<endl;
        cerr<<"P_dif\t"<<P_dif[t0]<<"\t"<<P_dif[t1]<<"\t"<<P_dif[ttau]<<endl;*/
//         cerr<<ftauArray<<endl;
//         cerr<<fArray<<endl;
//     }
}


// -----  Summary and intervention functions  -----

void MosqTransmission::uninfectVectors() {
    O_v.set_all( 0.0 );
    S_v.set_all( 0.0 );
    P_dif.set_all( 0.0 );
}

double sum1( const vecDay<double>& arr, SimTime end, SimTime N_v_length ){
    double val = 0.0;
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    for( SimTime d1 = end - SimTime::oneTS(); d1 < end; d1 += SimTime::oneDay() ){
        val += arr[mod_nn(d1, N_v_length)];
    }
    return val / SimTime::oneTS().inDays();
}
double sum2( const vecDay2D<double>& arr, SimTime end, SimTime N_v_length ){
    double val = 0.0;
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    for( SimTime d1 = end - SimTime::oneTS(); d1 < end; d1 += SimTime::oneDay() ){
        SimTime i1 = mod_nn(d1, N_v_length);
        for( size_t g = 0; g < Genotypes::N(); ++g ){
            val += arr.at(i1, g);
        }
    }
    return val / SimTime::oneTS().inDays();
}
double sum3( const vecDay2D<double>& arr, size_t g, SimTime end, SimTime N_v_length ){
    double val = 0.0;
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    for( SimTime d1 = end - SimTime::oneTS(); d1 < end; d1 += SimTime::oneDay() ){
        val += arr.at(mod_nn(d1, N_v_length), g);
    }
    return val / SimTime::oneTS().inDays();
}
double MosqTransmission::getLastVecStat( VecStat vs )const{
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    // One plus last, plus (0 mod N_v_length) to avoid negatives:
    SimTime end = sim::now() + SimTime::oneDay() + N_v_length;
    switch( vs ){
        case PA: return sum1(P_A, end, N_v_length);
        case PDF: return sum1(P_df, end, N_v_length);
        case PDIF: return sum2(P_dif, end, N_v_length);
        case NV: return sum1(N_v, end, N_v_length);
        case OV: return sum2(O_v, end, N_v_length);
        case SV: return sum2(S_v, end, N_v_length);
        default: throw SWITCH_DEFAULT_EXCEPTION;
    }
}
void MosqTransmission::summarize( size_t species )const{
    // Last time step ended at sim::now(). Values are stored per day, and for
    // the last time step values at sim::now() and four previos were set.
    // One plus last, plus (0 mod N_v_length) to avoid negatives:
    SimTime end = sim::now() + SimTime::oneDay() + N_v_length;
    mon::reportStatMSF( mon::MVF_LAST_NV0, species, getLastN_v0() );
    mon::reportStatMSF( mon::MVF_LAST_NV, species, sum1(N_v, end, N_v_length) );
    for( size_t g = 0; g < Genotypes::N(); ++g ){
        mon::reportStatMSGF( mon::MVF_LAST_OV, species, g, sum3(O_v, g, end, N_v_length) );
        mon::reportStatMSGF( mon::MVF_LAST_SV, species, g, sum3(S_v, g, end, N_v_length) );
    }
}

}
}
}
