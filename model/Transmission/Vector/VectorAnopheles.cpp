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

#include "Transmission/Vector/VectorAnopheles.h"
#include "Transmission/Vector/Nv0DelayFitting.h"
#include "Transmission/PerHostTransmission.h"
#include "Transmission/TransmissionModel.h"
#include "Host/Human.h"

#include "inputData.h"
#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include <fstream>
#include <ctime>

namespace OM {
namespace Transmission {
using namespace OM::util;


// -----  Initialisation of model, done before human warmup  ------

string VectorAnopheles::initialise (
    const scnXml::AnophelesParams& anoph,
    vector<double>& initialisationEIR,
    map<string, double>& nonHumanHostPopulations,
    int populationSize
)
{
    // -----  Set model variables  -----

    const scnXml::Mosq& mosq = anoph.getMosq();

    mosqRestDuration = mosq.getMosqRestDuration().getValue();
    EIPDuration = mosq.getExtrinsicIncubationPeriod().getValue();

    mosqSeekingDuration = mosq.getMosqSeekingDuration().getValue();

    probMosqSurvivalOvipositing = mosq.getMosqProbOvipositing().getValue();

    humanBase = mosq;

    minInfectedThreshold = mosq.getMinInfectedThreshold();

    if (1 > mosqRestDuration || mosqRestDuration > EIPDuration) {
        throw util::xml_scenario_error ("Code expects EIPDuration >= mosqRestDuration >= 1");
    }
    N_v_length = EIPDuration + mosqRestDuration;

    
    initAvailability( anoph, nonHumanHostPopulations, populationSize );

    initEIR( anoph, initialisationEIR );
    
    return anoph.getMosquito();
}

void VectorAnopheles::scaleEIR( double factor ) {
    FSCoeffic[0] += log( factor );
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

void VectorAnopheles::initAvailability(
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

double VectorAnopheles::calcEntoAvailability(double N_i, double P_A, double P_Ai)
{
    return (1.0 / N_i)
           * (P_Ai / (1.0-P_A))
           * (-log(P_A) / mosqSeekingDuration);
}

void VectorAnopheles::initEIR(
    const scnXml::AnophelesParams& anoph,
    vector<double>& initialisationEIR
){
    const scnXml::Seasonality& seasonality = anoph.getSeasonality();
    if( seasonality.getInput() != "EIR" ){
        throw util::xml_scenario_error("entomology.anopheles.seasonality.input: must be EIR (for now)");
        //TODO
    }
    // EIR for this species, with index 0 refering to value over first interval
    vector<double> speciesEIR (TimeStep::fromYears(1).inDays());

    if ( seasonality.getFourierSeries().present() ) {
        const scnXml::FourierSeries& seasFC = seasonality.getFourierSeries().get();
        const scnXml::FourierSeries::CoefficSequence& fsCoeffic = seasFC.getCoeffic();
        
        FSCoeffic.reserve (2*fsCoeffic.size() + 1);

        FSCoeffic.push_back( 0.0 );     // value doesn't matter; EIR will be scaled
        for( scnXml::FourierSeries::CoefficConstIterator it=fsCoeffic.begin(); it!=fsCoeffic.end(); ++it ){
            FSCoeffic.push_back( it->getA() );
            FSCoeffic.push_back( it->getB() );
        }
        // According to spec, EIR for first day of year (rather than EIR at the
        // exact start of the year) is generated with t=0 in Fourier series.
        EIRRotateAngle = seasFC.getEIRRotateAngle();
    } else if( seasonality.getMonthlyValues().present() ) {
        const scnXml::MonthlyValues& seasM = seasonality.getMonthlyValues().get();
        if( seasM.getSmoothing() != "fourier" ){
            throw util::xml_scenario_error("entomology.anopheles.seasonality.monthlyValues.smoothing: only fourier supported at the moment");
            //TODO: should be easy to add no smoothing
        }

        const size_t N_m = 12;
        const scnXml::MonthlyValues::ValueSequence seq = seasM.getValue();
        assert( seq.size() == N_m );    // enforced by schema
        double months[N_m];
        double sum = 0.0;
        for ( size_t i = 0; i < N_m; ++i ) {
            months[i] = seq[i];
            sum += months[i];
        }
        // arbitrary minimum we allow (cannot have zeros since we take the logarithm)
        double min = sum/1000.0;
        for ( size_t i = 0; i < N_m; ++i ) {
            if ( months[i] < min )
                months[i] = min;
        }

        const double PI = 3.14159265;
        const double w = 2.0 * PI / N_m;
        FSCoeffic.assign( 5, 0.0 );

        // Note: we use our values as the left-hand-side of our regions
        for ( size_t i = 0; i < N_m; ++i ) {
            double val = log( months[i] );
            FSCoeffic[0] += val;
            FSCoeffic[1] += val * cos( w*i );
            FSCoeffic[2] += val * sin( w*i );
            FSCoeffic[3] += val * cos( 2.0*w*i );
            FSCoeffic[4] += val * sin( 2.0*w*i );
        }
        FSCoeffic[0] /=N_m;
        FSCoeffic[1] *= 2.0 / N_m;
        FSCoeffic[2] *= 2.0 / N_m;
        FSCoeffic[3] *= 2.0 / N_m;
        FSCoeffic[4] *= 2.0 / N_m;

        // The above places the value for the first month at angle 0, so
        // effectively the first month starts at angle -2*pi/24 radians.
        // The value for the first day of the year should start 2*pi/(365*2)
        // radians later, so adjust EIRRotateAngle to compensate.
        EIRRotateAngle = M_PI * ( 1.0/12.0 - 1.0/365.0 );
    } else {
        assert( seasonality.getDailyValues().present() );      // XML loading code should enforce this
        throw util::xml_scenario_error("entomology.anopheles.seasonality.dailyValues: not supported yet");
        //TODO
    }
    
    if( !seasonality.getAnnualEIR().present() ){
        //TODO: work out when this is not required and implement code
        throw util::xml_scenario_error("entomology.anopheles.seasonality.annualEIR is required at the moment");
    }
    double targetEIR = seasonality.getAnnualEIR().get();
    
    // Now we rescale to get an EIR of targetEIR.
    // Calculate current sum as is usually done.
    calcFourierEIR (speciesEIR, FSCoeffic, EIRRotateAngle);
    // And scale:
    FSCoeffic[0] += log( targetEIR / vectors::sum( speciesEIR ) );

    // Calculate forced EIR for pre-intervention phase from FSCoeffic:
    calcFourierEIR (speciesEIR, FSCoeffic, EIRRotateAngle);

    // Add to the TransmissionModel's EIR, used for the initalization phase:
    for (int i = 0; i < TimeStep::fromYears(1).inDays(); ++i) {
        // index 1 of initialisationEIR corresponds to first period of year
        initialisationEIR[(1 + i / TimeStep::interval) % TimeStep::stepsPerYear] += speciesEIR[i];
    }

    if ( util::CommandLine::option( util::CommandLine::PRINT_ANNUAL_EIR ) ) {
        cout << "Annual EIR for "<<anoph.getMosquito()
             << ": "<<vectors::sum( speciesEIR )<<endl;
    }

    // Set other data used for mosqEmergeRate calculation:
    FSRotateAngle = EIRRotateAngle - (EIPDuration+10)/365.*2.*M_PI;       // usually around 20 days; no real analysis for effect of changing EIPDuration or mosqRestDuration
    initNvFromSv = 1.0 / anoph.getPropInfectious();
    initNv0FromSv = initNvFromSv * anoph.getPropInfected();       // temporarily use of initNv0FromSv

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

    N_v  .resize (N_v_length);
    O_v  .resize (N_v_length);
    S_v  .resize (N_v_length);
    P_A  .resize (N_v_length);
    P_df .resize (N_v_length);
    P_dif.resize (N_v_length);

    quinquennialS_v.assign (TimeStep::fromYears(5).inDays(), 0.0);
    forcedS_v.resize (TimeStep::fromYears(1).inDays());
    mosqEmergeRate.resize (TimeStep::fromYears(1).inDays()); // Only needs to be done here if loading from checkpoint
}


// -----  Initialisation of model which is done after running the human warmup  -----

void VectorAnopheles::setupNv0 (size_t sIndex,
                                const std::list<Host::Human>& population,
                                int populationSize,
                                double invMeanPopAvail) {
    // -----  N_v0, N_v, O_v, S_v  -----
    //BEGIN P_A, P_Ai, P_df, P_dif
    // rate at which mosquitoes find hosts or die (i.e. leave host-seeking state)
    double leaveSeekingStateRate = mosqSeekingDeathRate;

    // speciesEIR is average EIR per human over human population
    // that is, 1/populationSize * sum_{i in population} (P_Ai * P_B_i)
    // let sumPFindBite be sum_{i in population} (P_Ai * P_B_i):
    double sumPFindBite = 0.0;

    // NC's non-autonomous model provides two methods for calculating P_df and
    // P_dif; here we assume that P_E is constant.
    double intP_df = 0.0;

    for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
        const PerHostTransmission& host = h->perHostTransmission;
        double prod = host.entoAvailabilityFull (humanBase, sIndex, h->getAgeInYears(), invMeanPopAvail);
        leaveSeekingStateRate += prod;
        prod *= host.probMosqBiting(humanBase, sIndex);
        sumPFindBite += prod;
        intP_df += prod * host.probMosqResting(humanBase, sIndex);
    }

    for (vector<NHHParams>::const_iterator nhh = nonHumans.begin(); nhh != nonHumans.end(); ++nhh) {
        leaveSeekingStateRate += nhh->entoAvailability;
        intP_df += nhh->probCompleteCycle;
        // Note: in model, we do the same for intP_dif, except in this case it's
        // multiplied by infectiousness of host to mosquito which is zero.
    }

    // Probability of a mosquito not finding a host this day:
    double intP_A = exp(-leaveSeekingStateRate * mosqSeekingDuration);
    double P_Ai_base = (1.0 - intP_A) / leaveSeekingStateRate;
    sumPFindBite *= P_Ai_base;
    intP_df  *= P_Ai_base * probMosqSurvivalOvipositing;
    //END P_A, P_Ai, P_df, P_dif


    double initOvFromSv = initNv0FromSv;  // temporarily use of initNv0FromSv
    initNv0FromSv = initNvFromSv * (1.0 - intP_A - intP_df);

    // same as multiplying resultant eir since calcFourierEIR takes exp(...)
    FSCoeffic[0] += log (populationSize / sumPFindBite);
    calcFourierEIR (forcedS_v, FSCoeffic, FSRotateAngle);

    // Crude estimate of mosqEmergeRate: (1 - P_A(t) - P_df(t)) / (T * ρ_S) * S_T(t)
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);

    // Initialize per-day variables; S_v, N_v and O_v are only estimated
    for (int t = 0; t < N_v_length; ++t) {
        P_A[t] = intP_A;
        P_df[t] = intP_df;
        P_dif[t] = 0.0;     // humans start off with no infectiousness.. so just wait
        S_v[t] = forcedS_v[t];	// assume N_v_length ≤ 365
        N_v[t] = S_v[t] * initNvFromSv;
        O_v[t] = S_v[t] * initOvFromSv;
    }

    // All set up to drive simulation from forcedS_v
}


bool VectorAnopheles::vectorInitIterate () {
    // Try to match S_v against its predicted value. Don't try with N_v or O_v
    // because the predictions will change - would be chasing a moving target!
    // EIR comes directly from S_v, so should fit after we're done.

    double factor = vectors::sum (forcedS_v)*5 / vectors::sum(quinquennialS_v);
    //cout << "Pre-calced Sv, dynamic Sv:\t"<<sumAnnualForcedS_v<<'\t'<<vectors::sum(annualS_v)<<endl;
    if (!(factor > 1e-6 && factor < 1e6)) {
        if ( vectors::sum(forcedS_v) == 0.0 ) {
            return false;   // no EIR desired: nothing to do
        }
        cerr << "Input S_v for this vector:\t"<<vectors::sum(forcedS_v)<<endl;
        cerr << "Simulated S_v:\t\t\t"<<vectors::sum(quinquennialS_v)/5.0<<endl;
        throw util::traced_exception ("factor out of bounds (likely a code error)",util::Error::VectorFitting);
    }

    //cout << "Vector iteration: adjusting with factor "<<factor<<endl;
    // Adjusting mosqEmergeRate is the important bit. The rest should just
    // bring things to a stable state quicker.
    initNv0FromSv *= factor;
    initNvFromSv *= factor;     //(not currently used)
    vectors::scale (mosqEmergeRate, factor);
    vectors::scale (N_v, factor);
    // What factor exactly these should be scaled by isn't obvious; in any case
    // they should reach stable values quickly.
    vectors::scale (O_v, factor);
    vectors::scale (S_v, factor);
    vectors::scale (quinquennialS_v, factor); // scale so we can fit rotation offset

    // average annual period of S_v over 5 years
    vector<double> avgAnnualS_v( TimeStep::fromYears(1).inDays(), 0.0 );
    for ( int i = 0; i < TimeStep::fromYears(5).inDays(); ++i ) {
        avgAnnualS_v[i % TimeStep::fromYears(1).inDays()] =
            quinquennialS_v[i] / 5.0;
    }

    // Once the amplitude is approximately correct, we try to find a
    // rotation offset.
    double rAngle = Nv0DelayFitting::fit<double> (EIRRotateAngle, FSCoeffic, avgAnnualS_v);
    //cout << "Vector iteration: rotating with angle (in radians): " << rAngle << endl;
    // annualS_v was already rotated by old value of FSRotateAngle, so increment:
    FSRotateAngle -= rAngle;
    calcFourierEIR (forcedS_v, FSCoeffic, FSRotateAngle);
    // We use the stored initXxFromYy calculated from the ideal population age-structure (at init).
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);

    const double LIMIT = 0.1;
    return (fabs(factor - 1.0) > LIMIT) ||
           (rAngle > LIMIT * 2*M_PI / TimeStep::stepsPerYear);
}


// Every TimeStep::interval days:
void VectorAnopheles::advancePeriod (const std::list<Host::Human>& population,
                                     int populationSize,
                                     size_t sIndex,
                                     bool isDynamic,
                                     double invMeanPopAvail) {
    if (TimeStep::simulation > larvicidingEndStep) {
        larvicidingEndStep = TimeStep::future;
        larvicidingIneffectiveness = 1.0;
    }


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


    //BEGIN P_A, P_Ai, P_df, P_dif
    // rate at which mosquitoes find hosts or die (i.e. leave host-seeking state

    double leaveSeekingStateRate = mosqSeekingDeathRate;

    // NC's non-autonomous model provides two methods for calculating P_df and
    // P_dif; here we assume that P_E is constant.
    double intP_df = 0.0;
    double intP_dif = 0.0;
    for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
        const PerHostTransmission& host = h->perHostTransmission;
        double prod = host.entoAvailabilityFull (humanBase, sIndex, h->getAgeInYears(), invMeanPopAvail);
        leaveSeekingStateRate += prod;
        prod *= host.probMosqBiting(humanBase, sIndex)
                * host.probMosqResting(humanBase, sIndex);
        intP_df += prod;
        intP_dif += prod * h->probTransmissionToMosquito();
    }

    for (vector<NHHParams>::const_iterator nhh = nonHumans.begin(); nhh != nonHumans.end(); ++nhh) {
        leaveSeekingStateRate += nhh->entoAvailability;
        intP_df += nhh->probCompleteCycle;
        // Note: in model, we do the same for intP_dif, except in this case it's
        // multiplied by infectiousness of host to mosquito which is zero.
    }

    // Probability of a mosquito not finding a host this day:
    double intP_A = exp(-leaveSeekingStateRate * mosqSeekingDuration);
    double P_Ai_base = (1.0 - intP_A) / leaveSeekingStateRate;

    intP_df  *= P_Ai_base * probMosqSurvivalOvipositing;
    intP_dif *= P_Ai_base * probMosqSurvivalOvipositing;

    // Summed per day:
    partialEIR = 0.0;

    // The code within the for loop needs to run per-day, wheras the main
    // simulation uses TimeStep::interval day (currently 5 day) time steps.
    // The transmission for time-step t depends on the state during days
    // (t×(I-1)+1) through (t×I) where I is TimeStep::interval.
    int firstDay = TimeStep::simulation.inDays() - TimeStep::interval + 1;
    for (size_t i = 0; i < (size_t)TimeStep::interval; ++i) {
        // Warning: with x<0, x%y can be negative (depending on compiler); avoid x<0.
        // We add N_v_length so that ((dMod - x) >= 0) for (x <= N_v_length).
        size_t dMod = i + firstDay + N_v_length;
        assert (dMod >= (size_t)N_v_length);
        // Indecies for today, yesterday and mosqRestDuration days back:
        size_t t    = dMod % N_v_length;
        size_t t1   = (dMod - 1) % N_v_length;
        size_t ttau = (dMod - mosqRestDuration) % N_v_length;
        // Day of year and of 5-year cycles. Note that emergence during day 1
        // comes from mosqEmergeRate[0], hence subtraction by 1.
        size_t dYear1 = (firstDay + i - 1) % TimeStep::fromYears(1).inDays();
        size_t d5Year = (firstDay + i) % TimeStep::fromYears(5).inDays();


        // These only need to be calculated once per timestep, but should be
        // present in each of the previous N_v_length - 1 positions of arrays.
        P_A[t] = intP_A;
        P_df[t] = intP_df;
        P_dif[t] = intP_dif;


        //TODO: replace emergence with new formula: rho_p * (num pupae one day from emerging)
        // second two lines don't change
        N_v[t] = mosqEmergeRate[dYear1] * larvicidingIneffectiveness
                 + P_A[t1]  * N_v[t1]
                 + P_df[ttau] * N_v[ttau];
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
            sum += P_dif[tsl] * P_df[ttau] * (N_v[tsl] - O_v[tsl]) * ftauArray[EIPDuration+l-mosqRestDuration];
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


        if ( isDynamic ) {
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


        quinquennialS_v[d5Year] = S_v[t];

        partialEIR += S_v[t] * P_Ai_base;
    }
}


void VectorAnopheles::intervLarviciding (const scnXml::LarvicidingInterventionDescription& elt) {
    larvicidingIneffectiveness = 1 - elt.getEffectiveness().getValue();
    larvicidingEndStep = TimeStep::simulation + TimeStep(elt.getDuration().getValue());
}
void VectorAnopheles::uninfectVectors() {
    O_v.assign( O_v.size(), 0.0 );
    S_v.assign( S_v.size(), 0.0 );
    P_dif.assign( P_dif.size(), 0.0 );
}

double VectorAnopheles::getLastN_v0 (){
    double timestep_N_v0 = 0.0;
    int firstDay = TimeStep::simulation.inDays() - TimeStep::interval + 1;
    for (size_t i = 0; i < (size_t)TimeStep::interval; ++i) {
        size_t dYear1 = (firstDay + i - 1) % TimeStep::fromYears(1).inDays();
        timestep_N_v0 += mosqEmergeRate[dYear1];
    }
    return timestep_N_v0;
}
double VectorAnopheles::getLastVecStat ( VecStat vs ){
    //Note: implementation isn't performance optimal but rather intended to
    //keep code size low and have no overhead if not used.
    vector<double> *array;
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
    int firstDay = TimeStep::simulation.inDays() - TimeStep::interval + 1;
    for (size_t i = 0; i < (size_t)TimeStep::interval; ++i) {
        size_t t = (i + firstDay) % N_v_length;
        val += (*array)[t];
    }
    return val;
}
void VectorAnopheles::summarize (const string speciesName, Monitoring::Survey& survey) {
    survey.set_Vector_Nv0 (speciesName, getLastN_v0()/TimeStep::interval);
    survey.set_Vector_Nv (speciesName, getLastVecStat(NV)/TimeStep::interval);
    survey.set_Vector_Ov (speciesName, getLastVecStat(OV)/TimeStep::interval);
    survey.set_Vector_Sv (speciesName, getLastVecStat(SV)/TimeStep::interval);
}


void VectorAnopheles::calcFourierEIR (vector<double>& tArray, vector<double>& FC, double rAngle) {
    if (FC.size() % 2 == 0)
        throw util::xml_scenario_error("The number of Fourier coefficents should be odd.");

    // Frequency
    double w = 2*M_PI / double(tArray.size());

    // Number of Fourier Modes.
    int Fn = (FC.size()-1)/2;

    // Calculate inverse discrete Fourier transform
    for (size_t t=0; t<tArray.size(); t++) {
        double temp = FC[0];
        double wt = w*t - rAngle;
        for (int n=1;n<=Fn;n++) {
            temp = temp + FC[2*n-1]*cos(n*wt) + FC[2*n]*sin(n*wt);
        }
        tArray[t] = exp(temp);
    }
}

}
}
