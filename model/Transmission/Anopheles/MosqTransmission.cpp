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

//TODO: trim includes
#include "Transmission/Anopheles/MosqTransmission.h"
#include "Transmission/Anopheles/FixedEmergence.h"
#include "Transmission/Anopheles/SimpleMPDEmergence.h"
#include "schema/entomology.h"
#include "util/vectors.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;


// -----  Initialisation of model, done before human warmup  ------

MosqTransmission::MosqTransmission() :
        mosqRestDuration(0),
        EIPDuration(0),
        N_v_length(0),
        minInfectedThreshold( std::numeric_limits< double >::quiet_NaN() ),     // requires config
        timestep_N_v0(0.0)
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
        emergence = shared_ptr<EmergenceModel>( new LCEmergence() );
        emergence->initLifeCycle( lcOpt.get() );
        */
    }else if (util::ModelOptions::option( util::VECTOR_SIMPLE_MPD_MODEL )){
        if (!simpleMPDOpt.present())
            throw util::xml_scenario_error(
                "VECTOR_SIMPLE_MPD_MODEL: requires <simpleMPD> element with "
                "model parameters for each anopheles species");
        emergence = shared_ptr<EmergenceModel>(new SimpleMPDEmergence(simpleMPDOpt.get()) );
    }else
        emergence = shared_ptr<EmergenceModel>( new FixedEmergence() );
    
    
    // -----  Set model variables  -----

    mosqRestDuration = mosq.getMosqRestDuration().getValue();
    EIPDuration = mosq.getExtrinsicIncubationPeriod().getValue();
    if (1 > mosqRestDuration || mosqRestDuration*2 >= EIPDuration) {
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
    fArray.resize(EIPDuration-mosqRestDuration+1);
    fArray[0] = 1.0;
    ftauArray.resize(EIPDuration);
    for (int i = 0; i < mosqRestDuration; ++i)
        ftauArray[i] = 0.0;
    ftauArray[mosqRestDuration] = 1.0;
}

void MosqTransmission::initIterateScale ( double factor ){
    vectors::scale (N_v, factor);
    // What factor exactly these should be scaled by isn't obvious; in any case
    // they should reach stable values quickly.
    vectors::scale (O_v, factor);
    vectors::scale (S_v, factor);
}

void MosqTransmission::initState ( double tsP_A, double tsP_df,
                                       double initNvFromSv, double initOvFromSv,
                                       const vector<double>& forcedS_v ){
    N_v  .resize (N_v_length);
    O_v  .resize (N_v_length);
    S_v  .resize (N_v_length);
    P_A  .resize (N_v_length);
    P_df .resize (N_v_length);
    P_dif.resize (N_v_length);
    
    // Initialize per-day variables; S_v, N_v and O_v are only estimated
    assert( N_v_length <= static_cast<int>(forcedS_v.size()) );
    for (int t = 0; t < N_v_length; ++t) {
        P_A[t] = tsP_A;
        P_df[t] = tsP_df;
        P_dif[t] = 0.0;     // humans start off with no infectiousness.. so just wait
        //TODO: offset often won't be correct (comment from life-cycle code: is it true?)
        S_v[t] = forcedS_v[t];
        N_v[t] = S_v[t] * initNvFromSv;
        O_v[t] = S_v[t] * initOvFromSv;
    }
}

double MosqTransmission::update( size_t d, double tsP_A, double tsP_df,
                                 double tsP_dif, bool isDynamic, bool printDebug ){
    // Warning: with x<0, x%y can be negative (depending on compiler); avoid x<0.
    // We add N_v_length so that ((dMod - x) >= 0) for (x <= N_v_length).
    size_t dMod = d + N_v_length;
    assert (dMod >= (size_t)N_v_length);
    // Indecies for today, yesterday and mosqRestDuration days back:
    size_t t    = dMod % N_v_length;
    size_t t1   = (dMod - 1) % N_v_length;
    size_t ttau = (dMod - mosqRestDuration) % N_v_length;
    // Day of year. Note that emergence during day 1
    // comes from mosqEmergeRate[0], hence subtraction by 1.
    size_t dYear1 = (d - 1) % TimeStep::DAYS_IN_YEAR;
    
    // These only need to be calculated once per timestep, but should be
    // present in each of the previous N_v_length - 1 positions of arrays.
    P_A[t] = tsP_A;
    P_df[t] = tsP_df;
    P_dif[t] = tsP_dif;
    
    
    double nOvipositing = P_df[ttau] * N_v[ttau];
    double newAdults = emergence->get( d, dYear1, nOvipositing );
    util::streamValidate( newAdults );
    
    // num seeking mosquitos is: new adults + those which didn't find a host
    // yesterday + those who found a host tau days ago and survived cycle:
    N_v[t] = newAdults
                + P_A[t1]  * N_v[t1]
                + nOvipositing;
    // similar for O_v, except new mosquitoes are those who were uninfected
    // tau days ago, started a feeding cycle then, survived and got infected:
    O_v[t] = P_dif[ttau] * (N_v[ttau] - O_v[ttau])
                + P_A[t1]  * O_v[t1]
                + P_df[ttau] * O_v[ttau];
    
    //BEGIN S_v
    // Set up array with n in 1..θ_s−1 for f_τ(dMod-n) (NDEMD eq. 1.7)
    size_t fProdEnd = 2*mosqRestDuration;
    for (size_t n = mosqRestDuration+1; n <= fProdEnd; ++n) {
        size_t tn = (dMod-n)%N_v_length;
        ftauArray[n] = ftauArray[n-1] * P_A[tn];
    }
    ftauArray[fProdEnd] += P_df[(dMod-fProdEnd)%N_v_length];

    for (int n = fProdEnd+1; n < EIPDuration; ++n) {
        size_t tn = (dMod-n)%N_v_length;
        ftauArray[n] =
            P_df[tn] * ftauArray[n - mosqRestDuration]
            + P_A[tn] * ftauArray[n-1];
    }

    double sum = 0.0;
    size_t ts = dMod - EIPDuration;
    for (int l = 1; l < mosqRestDuration; ++l) {
        size_t tsl = (ts - l) % N_v_length;       // index dMod - theta_s - l
        sum += P_dif[tsl] * P_df[ttau] * (N_v[tsl] - O_v[tsl]) *
                ftauArray[EIPDuration+l-mosqRestDuration];
    }


    // Set up array with n in 1..θ_s−τ for f(dMod-n) (NDEMD eq. 1.6)
    for (int n = 1; n <= mosqRestDuration; ++n) {
        size_t tn = (dMod-n)%N_v_length;
        fArray[n] = fArray[n-1] * P_A[tn];
    }
    fArray[mosqRestDuration] += P_df[ttau];

    fProdEnd = EIPDuration-mosqRestDuration;
    for (size_t n = mosqRestDuration+1; n <= fProdEnd; ++n) {
        size_t tn = (dMod-n)%N_v_length;
        fArray[n] =
            P_df[tn] * fArray[n - mosqRestDuration]
            + P_A[tn] * fArray[n-1];
    }


    ts = ts % N_v_length;       // index dMod - theta_s
    S_v[t] = P_dif[ts] * fArray[EIPDuration-mosqRestDuration] * (N_v[ts] - O_v[ts])
                + sum
                + P_A[t1]*S_v[t1]
                + P_df[ttau]*S_v[ttau];


    if( isDynamic ){
        // We cut-off transmission when no more than X mosquitos are infected to
        // allow true elimination in simulations. Unfortunately, it may cause problems with
        // trying to simulate extremely low transmission, such as an R_0 case.
        if ( S_v[t] <= minInfectedThreshold ) { // infectious mosquito cut-off
            S_v[t] = 0.0;
            /* Note: could report; these reports often occur too frequently, however
            if( S_v[t] != 0.0 ){        // potentially reduce reporting
	cerr << TimeStep::simulation <<":\t S_v cut-off"<<endl;
            } */
        }
    }
    //END S_v
    
    //TODO: it should be possible to merge this call with emergence->get()
    // (S_v can be calculated before N_v).
    emergence->updateStats( d, tsP_dif, S_v[t] );
    
    timestep_N_v0 += newAdults;
    
    if( printDebug ){
        cerr<<"day "<<d<<":\temergence "<<newAdults<<",\tN_v "<<N_v[t]<<",\tS_v "<<S_v[t]<<endl;
        //cerr << "len: "<<N_v_length<<"\tdMod: "<<dMod<<"\tt: "<<t<<" "<<t1<<" "<<ttau<<"\tdYear1: "<<dYear1<<endl;
/*        cerr<<"P_A\t"<<P_A[t]<<"\t"<<P_A[t1]<<"\t"<<P_A[ttau]<<endl;
        cerr<<"P_df\t"<<P_df[t]<<"\t"<<P_df[t1]<<"\t"<<P_df[ttau]<<endl;
        cerr<<"P_dif\t"<<P_dif[t]<<"\t"<<P_dif[t1]<<"\t"<<P_dif[ttau]<<endl;*/
        cerr<<ftauArray<<endl;
        cerr<<fArray<<endl;
    }
    
    return S_v[t];
}


// -----  Summary and intervention functions  -----

void MosqTransmission::uninfectVectors() {
    O_v.assign( O_v.size(), 0.0 );
    S_v.assign( S_v.size(), 0.0 );
    P_dif.assign( P_dif.size(), 0.0 );
}

double MosqTransmission::getLastVecStat ( VecStat vs ) const{
    //Note: implementation isn't performance optimal but rather intended to
    //keep code size low and have no overhead if not used.
    const vector<double> *array = 0;
    switch( vs ){
        case PA: array = &P_A; break;
        case PDF: array = &P_df; break;
        case PDIF: array = &P_dif; break;
        case NV: array = &N_v; break;
        case OV: array = &O_v; break;
        case SV: array = &S_v; break;
        default: assert( false );
    }
    double val = 0.0;
    // add N_v_length to make sure LHS of % below is non-negative:
    int firstDay = TimeStep::simulation.inDays() - TimeStep::interval + 1 + N_v_length;
    for (size_t i = 0; i < (size_t)TimeStep::interval; ++i) {
        size_t t = (i + firstDay) % N_v_length;
        val += (*array)[t];
    }
    return val / TimeStep::interval;
}
void MosqTransmission::summarize (const string speciesName, Monitoring::Survey& survey) const{
    survey.set_Vector_Nv0 (speciesName, getLastN_v0());
    survey.set_Vector_Nv (speciesName, getLastVecStat(NV));
    survey.set_Vector_Ov (speciesName, getLastVecStat(OV));
    survey.set_Vector_Sv (speciesName, getLastVecStat(SV));
}

}
}
}
