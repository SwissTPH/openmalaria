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

#include "Global.h"
#include "Transmission/Anopheles/AnophelesModel.h"
#include "Transmission/PerHost.h"
#include "Population.h"
#include "WithinHost/Genotypes.h"
#include "util/vectors.h"
#include "util/errors.h"

#include <cmath>
#include <boost/format.hpp>

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;


// -----  Initialisation of model, done before human warmup  ------

string AnophelesModel::initialise (
    const scnXml::AnophelesParams& anoph,
    vector<double>& initialisationEIR,
//     map<string, double>& nonHumanHostPopulations,
    int populationSize
)
{
    // -----  Set model variables  -----

    const scnXml::Mosq& mosq = anoph.getMosq();

    mosqSeekingDuration = mosq.getMosqSeekingDuration().getValue();
    probMosqSurvivalOvipositing = mosq.getMosqProbOvipositing().getValue();
    humanBase = mosq;	// read human-specific parameters

    transmission.initialise( anoph.getLifeCycle(), anoph.getSimpleMPD(), anoph.getMosq() );
    
    // Uses anoph.getNonHumanHosts() and anoph.getMosq():
    initAvailability( anoph, /*nonHumanHostPopulations,*/ populationSize );
    
    // Uses anoph.getSeasonality() and three attributes:
    transmission.emergence->initEIR( anoph, initialisationEIR, transmission.getEIPDuration() );
    
    return anoph.getMosquito();
}


void AnophelesModel::initAvailability(
    const scnXml::AnophelesParams& anoph,
//     map<string, double>& nonHumanHostPopulations,
    int populationSize)
{
    // Set parameters. Notation is as in: Parameter Values for Transmission
    // Model, Chitnis et al. Sept 2010, equations (13), (14), (15).
    
    const scnXml::Mosq& mosq = anoph.getMosq();
    // A: Host seeking
    // Proportion of host-seeking parous mosquitoes (those who have laid eggs)
    // which laid eggs that day:
    const double A0 = mosq.getMosqLaidEggsSameDayProportion().getValue();
    // Probability that the mosquito survives the feeding cycle.
    // Note: Pf = M, the parous rate (prop mosqs which have laid eggs):
    const double Pf = mosq.getMosqSurvivalFeedingCycleProbability().getValue();
    const double humanBloodIndex = mosq.getMosqHumanBloodIndex().getValue();    // χ (chi)
    // Cycle probabilities, when biting a human:
    // B: Host encountered
    const double P_B1 = mosq.getMosqProbBiting().getMean();
    // C: Fed
    const double P_C1 = mosq.getMosqProbFindRestSite().getMean();
    // D: Resting
    const double P_D1 = mosq.getMosqProbResting().getMean();
    // E: Laying eggs (ovipositing)
    const double P_E1 = mosq.getMosqProbOvipositing().getValue();
    
    
    // -----  Calculate P_A, P_A1, P_A_n  -----
    // P_A is prob that a mosq is still host seeking after 1 day. It is also the
    // proportion of parous mosquitoes who have waited at least 1 day since
    // laying, thus 1 - P_A = A0.
    const double initP_A  = 1.0 - A0;
    // This is multiplied by the probability of encountering some type of host
    // on a given night is the total availability of this type of host (N_i * α_i).
    // Note: P_Ai = α_i * N_i / availFactor
    const double  availFactor = -log(initP_A) / (mosqSeekingDuration * (1.0 - initP_A));
    // Probability that a mosquito encounters a human on a given night
    double P_A1 = numeric_limits<double>::quiet_NaN();
    // Probability that a mosquito encounters any non human host on a given night
    // (confusingly labelled P_Ah in paper)
    double P_Ah = numeric_limits<double>::quiet_NaN();

    const scnXml::AnophelesParams::NonHumanHostsSequence& xmlSeqNNHs = anoph.getNonHumanHosts();
    if( xmlSeqNNHs.empty() ){
        // Number of non-human hosts: χ=1
        P_A1 = A0 * Pf / (P_B1 * P_C1 * P_D1 * P_E1);
        P_Ah = 0.0;
    }else{
        // Have non-human hosts: χ<1
        
        // let v = χ * P_D_1 * P_E_1; note that this is the average for humans
        const double v = humanBloodIndex * P_D1 * P_E1;
        const double chi1 = 1.0 - humanBloodIndex;    // 1 - χ
        
        double sum_xi = 0.0;    // sum rel. avail across NNHs; should be 1
        double sum_u = 0.0;     // sum u across NNHs where u = xi * P_B * P_C
        double sum_uvw = 0.0;   // sum u*(v+w) across NNHs where w = (1-χ)*P_D*P_E
        
        foreach( const scnXml::NonHumanHosts& xmlNNH, xmlSeqNNHs ){
            // availability population of hosts of this type relative to other non-human hosts:
            const double xi_i = xmlNNH.getMosqRelativeEntoAvailability().getValue();
            // cycle probabilities, when biting this type of host:
            const double P_B_i = xmlNNH.getMosqProbBiting().getValue();
            const double P_C_i = xmlNNH.getMosqProbFindRestSite().getValue();
            const double P_D_i = xmlNNH.getMosqProbResting().getValue();
            sum_xi += xi_i;
            const double u_i = xi_i * P_B_i * P_C_i;
            sum_u += u_i;
            const double w_i = chi1 * P_D_i * P_E1;     // note: we now assume P_E_i = P_E1
            sum_uvw += u_i * (v + w_i);
        }
        
        if ( !(sum_xi>0.9999 && sum_xi<1.0001) ){
            throw xml_scenario_error (
                "The sum of the relative entomological availability (ξ_i) "
                "across non-human hosts must be 1!"
            );
        }
        
        // equations (14), (15) of paper:
        P_A1 = (A0 * Pf * humanBloodIndex * sum_u) /
                (P_B1 * P_C1 * sum_uvw);
        P_Ah = (A0 * Pf * chi1) / sum_uvw;
    }
    
    
    // -----  Calculate availability rate of hosts (α_i) and non-human population data  -----
    humanBase.setEntoAvailability( (P_A1 / populationSize) * availFactor );
    
    nhh_avail = 0.0;
    nhh_sigma_df = 0.0;
    foreach( const scnXml::NonHumanHosts& xmlNNH, xmlSeqNNHs ){
//         map<string, double>::const_iterator pop = nonHumanHostPopulations.find(xmlNNH.getName());
//         if (pop == nonHumanHostPopulations.end()){
//             throw xml_scenario_error ((boost::format("There is no population size defined for "
//             "non-human host type \"%1%\"") %xmlNNH.getName()).str());
//         }
        
        // population size for non-human host category:
//         const double N_i = pop->second;
        // per-NNH parameters, as above:
        const double xi_i = xmlNNH.getMosqRelativeEntoAvailability().getValue();
        const double P_B_i = xmlNNH.getMosqProbBiting().getValue();
        const double P_C_i = xmlNNH.getMosqProbFindRestSite().getValue();
        const double P_D_i = xmlNNH.getMosqProbResting().getValue();
        const double P_Ahi = P_Ah * xi_i;       // probability of encountering this type of NNH on a given night
        const double avail_i = P_Ahi * availFactor; // N_i * α_i
        
        nhh_avail += avail_i;    // N * α
        nhh_sigma_df += avail_i * P_B_i * P_C_i * P_D_i;    // term in P_df series
        // Note: we would do the same for P_dif except that it's multiplied by
        // infectiousness of host to mosquito which is zero.
    }
    
    // ———  set mosqSeekingDeathRate  ———
    // Note: sum_k( P_Ak ) = P_A1 + P_Ah
    const double mu1 = (1.0 - (initP_A + P_A1 + P_Ah)) / (1.0 - initP_A);
    const double mu2 = -log(initP_A) / mosqSeekingDuration;
    mosqSeekingDeathRate = mu1 * mu2;
}


// -----  Initialisation of model which is done after creating initial humans  -----

void AnophelesModel::init2 (int nHumans, double meanPopAvail,
        double sum_avail, double sigma_f, double sigma_df)
{
    // -----  Calculate P_A, P_Ai, P_df based on pop age structure  -----
    
    // ν_A: rate at which mosquitoes find hosts or die (i.e. leave host-seeking state)
    double leaveRate = sum_avail + nhh_avail + mosqSeekingDeathRate;
    sigma_df += nhh_sigma_df;
    
    // Probability of a mosquito not finding a host this day:
    double initialP_A = exp(-leaveRate * mosqSeekingDuration);
    double availDivisor = (1.0 - initialP_A) / leaveRate;   // α_d
    // Input per-species EIR is the mean EIR experienced by a human adult.
    // We use sumPFindBite below to get required S_v.
    double sumPFindBite = sigma_f * availDivisor;
    double initialP_df  = sigma_df * availDivisor * probMosqSurvivalOvipositing;
    
    // -----  Calculate required S_v based on desired EIR  -----
    // Third parameter is a multiplication factor for S_v/EIR. First we multiply
    // input EIR by meanPopAvail to give us population average EIR instead of
    // adult average EIR, then we divide by (sumPFindBite/populationSize) to
    // get S_v.
    transmission.emergence->init2( initialP_A, initialP_df,
            nHumans * meanPopAvail / sumPFindBite, transmission );
    
    // All set up to drive simulation from forcedS_v
}

void AnophelesModel::initVectorInterv( const scnXml::VectorSpeciesIntervention& elt, size_t instance ){
    transmission.initVectorInterv( elt, instance );
    
    if( seekingDeathRateIntervs.size() <= instance )
        seekingDeathRateIntervs.resize( instance+1 );
    if( probDeathOvipositingIntervs.size() <= instance )
        probDeathOvipositingIntervs.resize( instance+1 );
    
    if( elt.getSeekingDeathRateIncrease().present() ){
        const scnXml::SeekingDeathRateIncrease& elt2 = elt.getSeekingDeathRateIncrease().get();
        if( elt2.getInitial() < -1.0 )
            throw util::xml_scenario_error( "seekingDeathRateIncrease intervention: initial effect must be ≥ -1" );
        seekingDeathRateIntervs[instance].set (elt2.getInitial(), elt2.getDecay(), "seekingDeathRateIncrease");
    }
    if( elt.getProbDeathOvipositing().present() ){
        const scnXml::ProbDeathOvipositing& elt2 = elt.getProbDeathOvipositing().get();
        if( elt2.getInitial() < 0.0 || elt2.getInitial() > 1.0 )
            throw util::xml_scenario_error( "probDeathOvipositing intrevention: initial effect must be in range [0,1]" );
        probDeathOvipositingIntervs[instance].set (elt2.getInitial(), elt2.getDecay(), "probDeathOvipositing");
    }
}
void AnophelesModel::initVectorTrap(const scnXml::Description1& desc, size_t instance){
    assert(trapParams.size() == instance);      // if triggered, this is a code error not XML
    TrapParams params;
    params.relAvail = desc.getRelativeAvailability().getValue();
    params.availDecay= DecayFunction::makeObject(desc.getDecayOfAvailability(), "decayOfAvailability");
    trapParams.push_back(params);
}

void AnophelesModel::deployVectorPopInterv (size_t instance){
    transmission.emergence->deployVectorPopInterv(instance);
    // do same as in above function (of EmergenceModel)
    assert( instance < seekingDeathRateIntervs.size() && instance < probDeathOvipositingIntervs.size() );
    seekingDeathRateIntervs[instance].deploy( sim::now() );
    probDeathOvipositingIntervs[instance].deploy( sim::now() );
}
void AnophelesModel::deployVectorTrap(size_t instance, double number, SimTime lifespan){
    assert(instance < trapParams.size());
    TrapData data;
    data.instance = instance;
    double adultAvail = humanBase.entoAvailability.mean();
    data.initialAvail = number * adultAvail * trapParams[instance].relAvail;
    data.availHet = trapParams[instance].availDecay->hetSample();
    data.deployTime = sim::now();
    data.expiry = sim::now() + lifespan;
    baitedTraps.push_back(data);
}

// Every SimTime::oneTS() days:
void AnophelesModel::advancePeriod (
        double sum_avail, double sigma_df, vector<double>& sigma_dif, bool isDynamic)
{
    transmission.emergence->update();
    
    /* Largely equations correspond to Nakul Chitnis's model in
      "A mathematic model for the dynamics of malaria in
      mosquitoes feeding on a heterogeneous host population" [MMDM]
    section 2, 3.5-3.6, plus extensions to a non-autonomous case from
      "Nonautonomous Difference Equations for Malaria Dynamics
                   in a Mosquito Population" [NDEMD]

    We calculate EIR over a time step (one or five days) as:
      sum_{for t over days in step} σ_i[t] * s_v[t]
      = sum_... (N_v[t] * P_Ai[t] * P_B_i[t])/(T*N_i[t]) * S_v[t]/N_v[t]
      = sum_... P_Ai[t] * P_B_i[t] * S_v[t]
    (since T == 1 and N_i[t] == 1 for all t).

      P_Ai[t] = (1 - P_A[t]) α_i[t] / sum_{h in hosts} α_h[t]
    (letting N_h[t] == 1 for all h,t). The only part of this varying per-host is
      α_i[t] = host.entoAvailability (index, human.getAgeInYears())
      Let availDivisor[t] = (1 - P_A[t]) / (sum_{h in hosts} α_h[t] + μ_vA).

    Note that although the model allows α_i and P_B_i to vary per-day, they only
    vary per time step of the main simulation. Hence:
      EIR = (sum_{t=...} S_v[t] * availDivisor[t]) * α_i * P_B_i

    Since S_v[t] * availDivisor[t] does not vary per individual, we calculate this
    per time step of the main simulation as partialEIR:
      partialEIR = (sum_{t=...} S_v[t] * availDivisor[t])

    Hence calculateEIR() only needs to do the following:
      EIR = partialEIR * α_i * P_B_i
    */


    // -----  Calculate P_A, P_Ai, P_df, P_dif based on human pop  -----
    
    // ν_A: rate at which mosquitoes find hosts or die (i.e. leave host-seeking state
    double leaveRate = mosqSeekingDeathRate;
    foreach( const util::SimpleDecayingValue& increase, seekingDeathRateIntervs ){
        leaveRate *= 1.0 + increase.current_value( sim::ts0() );
    }
    leaveRate += sum_avail;
    
    leaveRate += nhh_avail;
    sigma_df += nhh_sigma_df;
    
    for( list<TrapData>::iterator it = baitedTraps.begin(); it != baitedTraps.end(); ){
        if( sim::ts0() > it->expiry ){
            it = baitedTraps.erase(it);
            continue;
        }
        SimTime age = sim::ts0() - it->deployTime;
        double decayCoeff = trapParams[it->instance].availDecay->eval( age, it->availHet );
        leaveRate += it->initialAvail * decayCoeff;
        // sigma_df doesn't change: mosquitoes do not survive traps
        it++;
    }
    
    // Probability of a mosquito not finding a host this day:
    double tsP_A = exp(-leaveRate * mosqSeekingDuration);
    double availDivisor = (1.0 - tsP_A) / leaveRate;    // α_d
    
    // alphaE (α_E) is α_d * P_E, where P_E may be adjusted by interventions
    double alphaE = availDivisor * probMosqSurvivalOvipositing;
    foreach( const util::SimpleDecayingValue& pDeath, probDeathOvipositingIntervs ){
        alphaE *= 1.0 - pDeath.current_value( sim::ts0() );
    }
    double tsP_df  = sigma_df * alphaE;
    // from now, sigma_dif becomes P_dif (but we can't simply rename):
    vectors::scale( sigma_dif, alphaE );
    
    
    // Summed per day:
    partialEIR.assign( WithinHost::Genotypes::N(), 0.0 );
    
    transmission.resetTSStats();
    
    // The code within the for loop needs to run per-day, wheras the main
    // simulation uses one or five day time steps.
    const SimTime nextTS = sim::ts0() + SimTime::oneTS();
    for( SimTime d0 = sim::ts0(); d0 < nextTS; d0 += SimTime::oneDay() ){
        transmission.update( d0, tsP_A, tsP_df, sigma_dif, isDynamic, partialEIR, availDivisor );
    }
}

}
}
}
