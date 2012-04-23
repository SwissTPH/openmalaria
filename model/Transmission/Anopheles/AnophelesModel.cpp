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

#include "Transmission/Anopheles/AnophelesModel.h"
#include "Transmission/PerHost.h"
#include "Host/Human.h"

#include "util/errors.h"

#include <cmath>

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;


// -----  Initialisation of model, done before human warmup  ------

string AnophelesModel::initialise (
    const scnXml::AnophelesParams& anoph,
    vector<double>& initialisationEIR,
    map<string, double>& nonHumanHostPopulations,
    int populationSize
)
{
    // -----  Set model variables  -----

    const scnXml::Mosq& mosq = anoph.getMosq();

    mosqSeekingDuration = mosq.getMosqSeekingDuration().getValue();
    probMosqSurvivalOvipositing = mosq.getMosqProbOvipositing().getValue();
    humanBase = mosq;	// read human-specific parameters

    transmission.initialise( anoph );
    
    initAvailability( anoph, nonHumanHostPopulations, populationSize );
    
    transmission.emergence->initEIR( anoph, initialisationEIR, transmission.getEIPDuration() );
    
    return anoph.getMosquito();
}


/* Used by initAvailability to store temporary data. */
class AnophelesNonHumanParams {
public:
    void operator= (const scnXml::NonHumanHosts& nnh)
    {
        name = nnh.getName();
        relativeEntoAvailability = nnh.getMosqRelativeEntoAvailability().getValue();
        probMosqBiting = nnh.getMosqProbBiting().getValue();
        probMosqFindRestSite = nnh.getMosqProbFindRestSite().getValue();
        probMosqSurvivalResting = nnh.getMosqProbResting().getValue();
    }
    
    double prodBCD() const{	// P_B_i * P_C_i * P_D_i
        return probMosqBiting * probMosqFindRestSite * probMosqSurvivalResting;
    }
    
    string name;        // name of host type
    double relativeEntoAvailability;            // ξ_i
    double probMosqBiting;				// P_B_i
    double probMosqFindRestSite;		// P_C_i
    double probMosqSurvivalResting;	// P_D_i
};

void AnophelesModel::initAvailability(
    const scnXml::AnophelesParams& anoph,
    map<string, double>& nonHumanHostPopulations,
    int populationSize)
{
    // -----  Read XML data  -----
    
    const scnXml::AnophelesParams::NonHumanHostsSequence& otherHosts = anoph.getNonHumanHosts();
    vector<AnophelesNonHumanParams> nonHumanElts;
    nonHumanElts.resize (otherHosts.size());
    for (size_t i = 0; i < otherHosts.size(); ++i) {
        nonHumanElts[i] = otherHosts[i];
    }
    // read χ, P_B_1, P_C_1, P_D_1 and P_E_1; χ and P_E_1 are only needed to
    // calculate availability while the others are normally sampled.
    const scnXml::Mosq& mosq = anoph.getMosq();
    double mosqLaidEggsSameDayProp = mosq.getMosqLaidEggsSameDayProportion().getValue();
    double probMosqSurvivalFeedingCycle = mosq.getMosqSurvivalFeedingCycleProbability().getValue();
    double humanBloodIndex = mosq.getMosqHumanBloodIndex().getValue();
    double probBiting = mosq.getMosqProbBiting().getMean();
    double probFindRestSite = mosq.getMosqProbFindRestSite().getMean();
    double probResting = mosq.getMosqProbResting().getMean();
    double probOvipositing = mosq.getMosqProbOvipositing().getValue();
    
    
    // -----  Calculate P_A, P_A1, P_A_n  -----
    // Reference: Parameter Values for Transmission Model,
    // Chitnis et al. Sept 2010, equations (13), (14), (15).
    
    // Probability that a mosquito does not find a host and does not die in
    // one night of searching (P_A)
    double initP_A  = 1.0 - mosqLaidEggsSameDayProp;
    // Probability that a mosquito encounters a human on a given night
    double P_A1;
    // Probability that a mosquito encounters a non human host on a given night
    // (confusingly labelled P_Ah in paper)
    double P_Ah;

    if (nonHumanElts.empty()) {
        // i.e. χ=1
        
        // A_0 * P_f
        double pFedAndLaid = mosqLaidEggsSameDayProp * probMosqSurvivalFeedingCycle;
        // P_B_i * P_C_i * P_D_i * P_E_i (for average human)
        double pBiteRestOviposit = probBiting * probFindRestSite * probResting * probOvipositing;
        P_A1 = pFedAndLaid / pBiteRestOviposit;
        P_Ah = 0.0;
    }else{
        // i.e. χ<1
        
        // let v = χ * P_D_1 * P_E_1; note that this is the average for humans
        double v = humanBloodIndex * probResting * probOvipositing;
        // let chi1 = 1 - χ
        double chi1 = 1.0 - humanBloodIndex;
        
        // Sxi is the sum of ξ_i across non-human hosts i
        double Sxi = 0.0;
        // let u_i = ξ_i * P_B_i * P_C_i
        // Su is the sum of u_i across non-human hosts i
        double Su = 0.0;
        // let w_i = chi1 * P_D_i * P_E_i
        // Suvw is the sum of u_i*(v+w_i) across non-human hosts i
        double Suvw = 0.0;
        
        typedef vector<AnophelesNonHumanParams>::const_iterator ItType;
        for (ItType nhh = nonHumanElts.begin(); nhh != nonHumanElts.end(); ++nhh) {
            Sxi += nhh->relativeEntoAvailability;
            double u_i = nhh->relativeEntoAvailability * nhh->probMosqBiting * nhh->probMosqFindRestSite;
            Su += u_i;
            double w_i = chi1 * nhh->probMosqSurvivalResting * probOvipositing;
            Suvw += u_i * (v + w_i);
        }
        
        if ( !(Sxi>0.9999 && Sxi<1.0001) ){
            throw xml_scenario_error (
                "The sum of the relative entomological availability (ξ_i) "
                "across non-human hosts must be 1!"
            );
        }
        
        double A0Pf = mosqLaidEggsSameDayProp * probMosqSurvivalFeedingCycle;
        // P_A1 = A_0 * P_f * χ * Su  ...
        P_A1 = (A0Pf * humanBloodIndex * Su) /
        // ...  over P_B_1 * P_C_1 * Suvw
                (probBiting * probFindRestSite * Suvw);
        // and this one's as written
        P_Ah = (A0Pf * chi1) / Suvw;
    }
    
    
    // -----  Calculate availability rate of hosts (α_i) and non-human population data  -----
    humanBase.setEntoAvailability(
        calcEntoAvailability(populationSize, initP_A, P_A1));
    
    nonHumans.resize(nonHumanElts.size());
    for( size_t i=0; i<nonHumanElts.size(); ++i ){
        map<string, double>::const_iterator pop = nonHumanHostPopulations.find(nonHumanElts[i].name);
        if (pop == nonHumanHostPopulations.end()){
            throw xml_scenario_error ("There is no population size defined for "
            "at least one non human host type, please check the scenario file.");
        }
        
        double nonHumanPopulationSize = pop->second;
        double entoAvailability = calcEntoAvailability(nonHumanPopulationSize,
                                                       initP_A,
                                                       P_Ah*nonHumanElts[i].relativeEntoAvailability);
        
        nonHumans[i].entoAvailability = entoAvailability;
        nonHumans[i].probCompleteCycle = entoAvailability * nonHumanElts[i].prodBCD();
    }
    
    // -----  Calculate death rate while seeking (µ_vA)  -----
    // since sum_i(ξ_i)=1, sum_k(P_A_k)=P_A1+P_Ah
    double mu1 = (1.0-initP_A-P_A1-P_Ah) / (1.-initP_A);
    double mu2 = -log(initP_A) / mosqSeekingDuration;
    mosqSeekingDeathRate = mu1 * mu2;
}

double AnophelesModel::calcEntoAvailability(double N_i, double P_A, double P_Ai)
{
    return (1.0 / N_i)
           * (P_Ai / (1.0-P_A))
           * (-log(P_A) / mosqSeekingDuration);
}


// -----  Initialisation of model which is done after creating initial humans  -----

void AnophelesModel::init2 (size_t sIndex,
                                const std::list<Host::Human>& population,
                                int populationSize,
                                double meanPopAvail)
{
    // -----  Calculate P_A, P_Ai, P_df based on pop age structure  -----
    
    // rate at which mosquitoes find hosts or die (i.e. leave host-seeking state)
    double leaveSeekingStateRate = mosqSeekingDeathRate;

    // Input per-species EIR is the mean EIR experienced by a human adult.
    // We use sumPFindBite below to get required S_v.
    // Let sumPFindBite be sum_{i in population} (P_Ai * P_B_i):
    double sumPFindBite = 0.0;

    // NC's non-autonomous model provides two methods for calculating P_df and
    // P_dif; here we assume that P_E is constant.
    initialP_df = 0.0;

    for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
        const OM::Transmission::PerHost& host = h->perHostTransmission;
        double prod = host.entoAvailabilityFull (humanBase, sIndex, h->getAgeInYears());
        leaveSeekingStateRate += prod;
        prod *= host.probMosqBiting(humanBase, sIndex);
        sumPFindBite += prod;
        initialP_df += prod * host.probMosqResting(humanBase, sIndex);
    }

    for (vector<NHHParams>::const_iterator nhh = nonHumans.begin(); nhh != nonHumans.end(); ++nhh) {
        leaveSeekingStateRate += nhh->entoAvailability;
        initialP_df += nhh->probCompleteCycle;
        // Note: in model, we do the same for tsP_dif, except in this case it's
        // multiplied by infectiousness of host to mosquito which is zero.
    }

    // Probability of a mosquito not finding a host this day:
    initialP_A = exp(-leaveSeekingStateRate * mosqSeekingDuration);
    double P_Ai_base = (1.0 - initialP_A) / leaveSeekingStateRate;
    sumPFindBite *= P_Ai_base;
    initialP_df  *= P_Ai_base * probMosqSurvivalOvipositing;
    
    
    // -----  Calculate required S_v based on desired EIR  -----
    // Last parameter is a multiplication factor for S_v/EIR. First we multiply
    // input EIR by meanPopAvail to give us population average EIR instead of
    // adult average EIR, then we divide by (sumPFindBite/populationSize) to
    // get S_v.
    transmission.emergence->init2( initialP_A, initialP_df, populationSize * meanPopAvail / sumPFindBite, transmission );
    
    // All set up to drive simulation from forcedS_v
}


// Every TimeStep::interval days:
void AnophelesModel::advancePeriod (const std::list<Host::Human>& population,
                                     int populationSize,
                                     size_t sIndex,
                                     bool isDynamic) {
    transmission.emergence->update();
    
    /* Largely equations correspond to Nakul Chitnis's model in
      "A mathematic model for the dynamics of malaria in
      mosquitoes feeding on a heterogeneous host population" [MMDM]
    section 2, 3.5-3.6, plus extensions to a non-autonomous case from
      "Nonautonomous Difference Equations for Malaria Dynamics
                   in a Mosquito Population" [NDEMD]

    We calculate EIR over a 5-day TimeStep::interval as:
      sum_{for t over days} σ_i[t] * s_v[t]
      = sum_... (N_v[t] * P_Ai[t] * P_B_i[t])/(T*N_i[t]) * S_v[t]/N_v[t]
      = sum_... P_Ai[t] * P_B_i[t] * S_v[t]
    (since T == 1 and N_i[t] == 1 for all t).

      P_Ai[t] = (1 - P_A[t]) α_i[t] / sum_{h in hosts} α_h[t]
    (letting N_h[t] == 1 for all h,t). The only part of this varying per-host is
      α_i[t] = host.entoAvailability (index, h->getAgeInYears())
      Let P_Ai_base[t] = (1 - P_A[t]) / (sum_{h in hosts} α_h[t] + μ_vA).

    Note that although the model allows α_i and P_B_i to vary per-day, they only
    vary per TimeStep::interval of the main simulation. Hence:
      EIR = (sum_{t=...} S_v[t] * P_Ai_base[t]) * α_i * P_B_i

    Since S_v[t] * P_Ai_base[t] does not vary per individual, we calculate this
    per TimeStep::interval of the main simulation as partialEIR:
      partialEIR = (sum_{t=...} S_v[t] * P_Ai_base[t])

    Hence calculateEIR() only needs to do the following:
      EIR = partialEIR * α_i * P_B_i
    */


    // -----  Calculate P_A, P_Ai, P_df, P_dif based on human pop  -----
    
    // rate at which mosquitoes find hosts or die (i.e. leave host-seeking state
    double leaveSeekingStateRate = mosqSeekingDeathRate;

    // NC's non-autonomous model provides two methods for calculating P_df and
    // P_dif; here we assume that P_E is constant.
    double tsP_df = 0.0;
    double tsP_dif = 0.0;
    for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
        const OM::Transmission::PerHost& host = h->perHostTransmission;
        double prod = host.entoAvailabilityFull (humanBase, sIndex, h->getAgeInYears());
        leaveSeekingStateRate += prod;
        prod *= host.probMosqBiting(humanBase, sIndex)
                * host.probMosqResting(humanBase, sIndex);
        tsP_df += prod;
        tsP_dif += prod * h->probTransmissionToMosquito();
    }

    for (vector<NHHParams>::const_iterator nhh = nonHumans.begin(); nhh != nonHumans.end(); ++nhh) {
        leaveSeekingStateRate += nhh->entoAvailability;
        tsP_df += nhh->probCompleteCycle;
        // Note: in model, we do the same for tsP_dif, except in this case it's
        // multiplied by infectiousness of host to mosquito which is zero.
    }

    // Probability of a mosquito not finding a host this day:
    double tsP_A = exp(-leaveSeekingStateRate * mosqSeekingDuration);
    double P_Ai_base = (1.0 - tsP_A) / leaveSeekingStateRate;

    tsP_df  *= P_Ai_base * probMosqSurvivalOvipositing;
    tsP_dif *= P_Ai_base * probMosqSurvivalOvipositing;
    
    
    if( true ){//TODO: may use this condition soon...
        // Summed per day:
        partialEIR = 0.0;
        
        transmission.resetTSStats();
        
        // The code within the for loop needs to run per-day, wheras the main
        // simulation uses TimeStep::interval day (currently 5 day) time steps.
        // The transmission for time-step t depends on the state during days
        // (t×(I-1)+1) through (t×I) where I is TimeStep::interval.
        int firstDay = TimeStep::simulation.inDays() - TimeStep::interval + 1;
        for (size_t i = 0; i < (size_t)TimeStep::interval; ++i) {
            partialEIR += transmission.update( i + firstDay, tsP_A, tsP_df, tsP_dif, isDynamic, false ) * P_Ai_base;
        }
    }
}

}
}
}
