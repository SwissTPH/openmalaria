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
#include "Transmission/Anopheles/FixedEmergence.h"
#include "Transmission/Anopheles/SimpleMPDEmergence.h"
#include "Population.h"
#include "WithinHost/Genotypes.h"
#include "util/vectors.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"

#include <cmath>

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;
using WithinHost::Genotypes;

// -----  Initialisation of model, done before human warmup  ------

string AnophelesModel::initialise (size_t species, const scnXml::AnophelesParams& anoph, vector<double>& initialisationEIR, int populationSize)
{
    // -----  Set model variables  -----
    const scnXml::Mosq& mosq = anoph.getMosq();
    const scnXml::AnophelesParams::SimpleMPDOptional& simpleMPDOpt = anoph.getSimpleMPD();
    // const scnXml::AnophelesParams::LifeCycleOptional& lcOpt = anoph.getLifeCycle();

    mosqSeekingDuration = mosq.getMosqSeekingDuration().getValue();
    probMosqSurvivalOvipositing = mosq.getMosqProbOvipositing().getValue();

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

    // Uses anoph.getNonHumanHosts() and anoph.getMosq():
    initAvailability( species, anoph, /*nonHumanHostPopulations,*/ populationSize );
    
    // Uses anoph.getSeasonality() and three attributes:
    initEIR( anoph, initialisationEIR, getEIPDuration() );
    
    return anoph.getMosquito();
}


void AnophelesModel::initAvailability(
    size_t species,
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
        
        for( const scnXml::NonHumanHosts& xmlNNH : xmlSeqNNHs ){
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
    PerHostAnophParams::scaleEntoAvailability( species, (P_A1 / populationSize) * availFactor );
    
    nhh_avail = 0.0;
    nhh_sigma_df = 0.0;
    nhh_sigma_dff = 0.0;
    for( const scnXml::NonHumanHosts& xmlNNH : xmlSeqNNHs ){
        const double xi_i = xmlNNH.getMosqRelativeEntoAvailability().getValue();
        const double P_B_i = xmlNNH.getMosqProbBiting().getValue();
        const double P_C_i = xmlNNH.getMosqProbFindRestSite().getValue();
        const double P_D_i = xmlNNH.getMosqProbResting().getValue();
        const double rel_fecundity = (xmlNNH.getHostFecundityFactor().present() ? xmlNNH.getHostFecundityFactor().get().getValue() : 1.0);
        const double P_Ahi = P_Ah * xi_i;       // probability of encountering this type of NNH on a given night
        const double avail_i = P_Ahi * availFactor; // N_i * α_i
        
        nhh_avail += avail_i;    // N * α (nhh_avail is not used anymore)
        const double df = avail_i * P_B_i * P_C_i * P_D_i;    // term in P_df series
        nhh_sigma_df += df;
        nhh_sigma_dff += df * rel_fecundity;
        // Note: we would do the same for P_dif except that it's multiplied by
        // infectiousness of host to mosquito which is zero.

        string name = xmlNNH.getName();

        NHH nhh;
        nhh.avail_i = avail_i;
        nhh.P_B_I = P_B_i;
        nhh.P_C_I = P_C_i;
        nhh.P_D_I = P_D_i;
        nhh.rel_fecundity = rel_fecundity;
        nhh.expiry = SimTime::future();

        initNhh[name] = nhh;
    }
    
    // ———  set mosqSeekingDeathRate  ———
    // Note: sum_k( P_Ak ) = P_A1 + P_Ah
    const double mu1 = (1.0 - (initP_A + P_A1 + P_Ah)) / (1.0 - initP_A);
    const double mu2 = -log(initP_A) / mosqSeekingDuration;
    mosqSeekingDeathRate = mu1 * mu2;
}

/** Called by initialise function to init variables directly related to EIR
 * 
 * @param anoph Data from XML
 * @param initialisationEIR In/out parameter: TransmissionModel::initialisationEIR
 * @param EIPDuration parameter from MosqTransmission (used for an estimation)
 */
void AnophelesModel::initEIR(const scnXml::AnophelesParams& anoph, vector<double>& initialisationEIR, SimTime EIPDuration)
{
    const scnXml::Seasonality& seasonality = anoph.getSeasonality();
    if ( seasonality.getInput() != "EIR" ) {
        throw util::xml_scenario_error("entomology.anopheles.seasonality.input: must be EIR (for now)");
        //TODO
    }
    // EIR for this species, with index 0 refering to value over first interval
    vecDay<double> speciesEIR (SimTime::oneYear());

    if ( seasonality.getFourierSeries().present() ) {
        const scnXml::FourierSeries& seasFC = seasonality.getFourierSeries().get();
        const scnXml::FourierSeries::CoefficSequence& fsCoeffic = seasFC.getCoeffic();

        FSCoeffic.reserve (2*fsCoeffic.size() + 1);

        FSCoeffic.push_back( 0.0 );     // value doesn't matter; EIR will be scaled
        for( auto it=fsCoeffic.begin(); it!=fsCoeffic.end(); ++it ) {
            FSCoeffic.push_back( it->getA() );
            FSCoeffic.push_back( it->getB() );
        }
        // According to spec, EIR for first day of year (rather than EIR at the
        // exact start of the year) is generated with t=0 in Fourier series.
        EIRRotateAngle = seasFC.getEIRRotateAngle();
    } else if ( seasonality.getMonthlyValues().present() ) {
        const scnXml::MonthlyValues& seasM = seasonality.getMonthlyValues().get();
        if ( seasM.getSmoothing() != "fourier" ) {
            throw util::xml_scenario_error("entomology.anopheles.seasonality.monthlyValues.smoothing: only fourier supported at the moment");
            //TODO: should be easy to add no smoothing
        }

        const size_t N_m = 12;
        const scnXml::MonthlyValues::ValueSequence seq = seasM.getValue();
        assert( seq.size() == N_m );    // enforced by schema
        vector<double> months(N_m);
        double sum = 0.0;
        for( size_t i = 0; i < N_m; ++i ) {
            months[i] = seq[i];
            sum += months[i];
        }
        // arbitrary minimum we allow (cannot have zeros since we take the logarithm)
        double min = sum/1000.0;
        for( size_t i = 0; i < N_m; ++i ) {
            if ( months[i] < min )
                months[i] = min;
        }

        FSCoeffic.assign( 5, 0.0 );
        //TODO: determine whether to use Fourier Series Coefficient method or
        // Discrete Fourier Transform. Former is designed for integrals (which
        // roughly what we have?), the latter for discrete values. DFT doesn't
        // work well when number of intervals changes?
        vectors::logFourierCoefficients(months, FSCoeffic);

        // The above places the value for the first month at angle 0, so
        // effectively the first month starts at angle -2*pi/24 radians.
        // The value for the first day of the year should start 2*pi/(365*2)
        // radians later, so adjust EIRRotateAngle to compensate.

        // Change this to 0 later?
        EIRRotateAngle = M_PI * ( 1.0/12.0 - 1.0/365.0 );
    } else {
        assert( seasonality.getDailyValues().present() );      // XML loading code should enforce this
        throw util::xml_scenario_error("entomology.anopheles.seasonality.dailyValues: not supported yet");
        //TODO: can we encode as a Fourier series without smoothing?
    }

    if ( !seasonality.getAnnualEIR().present() ) {
        //TODO: work out when this is not required and implement code
        throw util::xml_scenario_error("entomology.anopheles.seasonality.annualEIR is required at the moment");
    }
    double targetEIR = seasonality.getAnnualEIR().get();
    
    // Now we rescale to get an EIR of targetEIR.
    // Calculate current sum as is usually done.
    vectors::expIDFT (speciesEIR, FSCoeffic, EIRRotateAngle);
    // And scale (also acts as a unit conversion):
    FSCoeffic[0] += log( targetEIR / vectors::sum( speciesEIR ) );
    
    // Calculate forced EIR for pre-intervention phase from FSCoeffic:
    vectors::expIDFT (speciesEIR, FSCoeffic, EIRRotateAngle);
    
    // Add to the TransmissionModel's EIR, used for the initalization phase.
    // Note: sum stays the same, units changes to per-time-step.
    for( SimTime i = SimTime::zero(); i < SimTime::oneYear(); i += SimTime::oneDay() ){
        // index 0 of initialisationEIR corresponds to first period of year
        initialisationEIR[mod_nn(i.inSteps(), sim::stepsPerYear())] += speciesEIR[i];
    }
    
    if ( util::CommandLine::option( util::CommandLine::PRINT_ANNUAL_EIR ) ) {
        cout << "Annual EIR for "<<anoph.getMosquito()
             << ": "<<vectors::sum( speciesEIR )<<endl;
    }

    // Set other data used for mosqEmergeRate calculation:
    FSRotateAngle = EIRRotateAngle - (EIPDuration.inDays()+10)/365.*2.*M_PI;       // usually around 20 days; no real analysis for effect of changing EIPDuration or mosqRestDuration
    initNvFromSv = 1.0 / anoph.getPropInfectious();
    initOvFromSv = initNvFromSv * anoph.getPropInfected();
}

// -----  Initialisation of model which is done after creating initial humans  -----

void AnophelesModel::init2 (int nHumans, double meanPopAvail,
        double sum_avail, double sigma_f, double sigma_df, double sigma_dff)
{
    // -----  Calculate P_A, P_Ai, P_df based on pop age structure  -----
    
    // ν_A: rate at which mosquitoes find hosts or die (i.e. leave host-seeking state)
    double leaveRate = sum_avail + nhh_avail + mosqSeekingDeathRate;
    sigma_df += nhh_sigma_df;
    sigma_dff += nhh_sigma_dff;
    
    // Probability of a mosquito not finding a host this day:
    double tsP_A = exp(-leaveRate * mosqSeekingDuration);
    double availDivisor = (1.0 - tsP_A) / leaveRate;   // α_d
    // Input per-species EIR is the mean EIR experienced by a human adult.
    // We use sumPFindBite below to get required S_v.
    double sumPFindBite = sigma_f * availDivisor;
    double tsP_df  = sigma_df * availDivisor * probMosqSurvivalOvipositing;
    double tsP_dff  = sigma_dff * availDivisor * probMosqSurvivalOvipositing;
    
    double initialP_Amu = (1-initialP_A) * mosqSeekingDeathRate/(mosqSeekingDeathRate + sum_avail + nhh_avail);
    double initialP_A1 = (1-initialP_A) * sum_avail/(mosqSeekingDeathRate + sum_avail + nhh_avail);
    double initialP_Ah = 0.0;
    for( auto it = initNhh.begin(); it != initNhh.end(); ++it){
        initialP_Ah += (1-initialP_A) * it->second.avail_i / (mosqSeekingDeathRate + sum_avail + nhh_avail);
    }

    // -----  Calculate required S_v based on desired EIR  -----
    initNv0FromSv = initNvFromSv * (1.0 - tsP_A - tsP_df);

    // We scale FSCoeffic to give us S_v instead of EIR.
    // Log-values: adding log is same as exponentiating, multiplying and taking
    // the log again.

    double EIRtoS_v = nHumans * meanPopAvail / sumPFindBite;

    FSCoeffic[0] += log( EIRtoS_v );
    vectors::expIDFT (forcedS_v, FSCoeffic, EIRRotateAngle);
    
    N_v  .assign (N_v_length, numeric_limits<double>::quiet_NaN());
    O_v  .assign (N_v_length, Genotypes::N(), numeric_limits<double>::quiet_NaN());
    S_v  .assign (N_v_length, Genotypes::N(), numeric_limits<double>::quiet_NaN());
    P_A  .assign (N_v_length, numeric_limits<double>::quiet_NaN());
    P_df .assign (N_v_length, numeric_limits<double>::quiet_NaN());
    P_dif.assign (N_v_length, Genotypes::N(), 0.0);// humans start off with no infectiousness.. so just wait
    P_dff.assign (N_v_length, numeric_limits<double>::quiet_NaN());
    P_Amu  .assign (N_v_length, numeric_limits<double>::quiet_NaN());
    P_A1  .assign (N_v_length, numeric_limits<double>::quiet_NaN());
    P_Ah  .assign (N_v_length, numeric_limits<double>::quiet_NaN());

    // Initialize per-day variables; S_v, N_v and O_v are only estimated
    assert( N_v_length <= forcedS_v.size() );
    for( SimTime t = SimTime::zero(); t < N_v_length; t += SimTime::oneDay() ){
        P_A[t] = tsP_A;
        P_Amu[t] = tsP_Amu;
        P_A1[t] = tsP_A1;
        P_Ah[t] = tsP_Ah;
        P_df[t] = tsP_df;
        P_dff[t] = tsP_dff;
        N_v[t] = forcedS_v[t] * initNvFromSv;
        for( size_t genotype = 0; genotype < Genotypes::N(); ++genotype ){
            S_v.at(t, genotype) = forcedS_v[t] * Genotypes::initialFreq(genotype);
            O_v.at(t,genotype) = S_v.at(t,genotype) * initOvFromSv;
        }
    }

    // Crude estimate of mosqEmergeRate: (1 - P_A(t) - P_df(t)) / (T * ρ_S) * S_T(t)
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);
    
    // All set up to drive simulation from forcedS_v

    scaleFactor = 1.0;
    shiftAngle = FSRotateAngle;
    scaled = false;
    rotated = false;

    emergence->init2(tsP_dff, initNvFromSv, forcedS_v, mosqEmergeRate, mosqRestDuration);
    // All set up to drive simulation from forcedS_v
}

bool AnophelesModel::initIterate ()
{
    vecDay<double> avgAnnualS_v( SimTime::oneYear(), 0.0 );
    for( SimTime i = SimTime::fromYearsI(4); i < SimTime::fromYearsI(5); i += SimTime::oneDay() ){
        avgAnnualS_v[mod_nn(i, SimTime::oneYear())] =
            quinquennialS_v[i];
    }

    double factor = vectors::sum(forcedS_v) / vectors::sum(avgAnnualS_v);
    
    //cout << "check: " << vectors::sum(forcedS_v) << " " << vectors::sum(avgAnnualS_v) << endl;
    //cout << "Pre-calced Sv, dynamic Sv:\t"<<sumAnnualForcedS_v<<'\t'<<vectors::sum(annualS_v)<<endl;
    if (!(factor > 1e-6 && factor < 1e6)) {
        if( factor > 1e6 && vectors::sum(quinquennialS_v) < 1e-3 ){
            throw util::base_exception("Simulated S_v is approx 0 (i.e.\
 mosquitoes are not infectious, before interventions). Simulator cannot handle this; perhaps\
 increase EIR or change the entomology model.", util::Error::VectorFitting);
        }
        if ( vectors::sum(forcedS_v) == 0.0 ) {
            return false;   // no EIR desired: nothing to do
        }
        cerr << "Input S_v for this vector:\t"<<vectors::sum(forcedS_v)<<endl;
        cerr << "Simulated S_v:\t\t\t"<<vectors::sum(quinquennialS_v)/5.0<<endl;
        throw TRACED_EXCEPTION ("factor out of bounds",util::Error::VectorFitting);
    }

    const double LIMIT = 0.1;

    if(fabs(factor - 1.0) > LIMIT)
    {
        scaled = false;
        double factorDiff = (scaleFactor * factor - scaleFactor) * 1.0;
        scaleFactor += factorDiff;
    }
    else
        scaled = true;

    double rAngle = findAngle(EIRRotateAngle, FSCoeffic, avgAnnualS_v);
    shiftAngle += rAngle;
    rotated = true;

    // cout << "EIRRotateAngle: " << EIRRotateAngle << " rAngle = " << rAngle << ", angle = " << shiftAngle << " scalefactor: " << scaleFactor << " , factor: " << factor << endl;

    // Compute forced_sv from the Fourrier Coeffs
    // shiftAngle rotate the vector to correct the offset between simulated and input EIR
    // shiftAngle is the offset between the 
    vectors::expIDFT(mosqEmergeRate, FSCoeffic, -shiftAngle);
    // Scale the vector according to initNv0FromSv to get the mosqEmergerate
    // scaleFactor scales the vector to correct the ratio between simulated and input EIR
    vectors::scale (mosqEmergeRate, scaleFactor * initNv0FromSv);

    //initNvFromSv *= scaleFactor;     //(not currently used)

    // What factor exactly these should be scaled by isn't obvious; in any case
    // they should reach stable values quickly.
    vectors::scale (N_v, factor);
    vectors::scale (O_v, factor);
    vectors::scale (S_v, factor);

    emergence->initIterate(factor, mosqEmergeRate);

    return !(scaled && rotated);
}
//@}

void AnophelesModel::initVectorInterv( const scnXml::VectorSpeciesIntervention& elt, size_t instance ){
    if( emergenceReduction.size() <= instance )
            emergenceReduction.resize( instance+1 );
        
        if( elt.getEmergenceReduction().present() ){
            const scnXml::EmergenceReduction& elt2 = elt.getEmergenceReduction().get();
            if( elt2.getInitial() > 1.0 )
                throw util::xml_scenario_error( "emergenceReduction intervention: initial effect must be ≤ 1" );
            emergenceReduction[instance].set (elt2.getInitial(), elt2.getDecay(), "emergenceReduction");
    }
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
    trapParams.push_back(move(params));
}
void AnophelesModel::initNonHumanHostsInterv(const scnXml::NonHumanHostsSpeciesIntervention& elt, const scnXml::DecayFunction& decay, size_t instance, string name ){
    if( elt.getReduceAvailability().present() ){
        const scnXml::ReduceAvailability& elt2 = elt.getReduceAvailability().get();
        if( elt2.getInitial() < -1.0 )
            throw util::xml_scenario_error( "reduceAvailability intervention: initial effect must be ≥ -1" );
        reduceNHHAvailability[name].resize(instance+1);
    if( reduceP_B_I[name].size() <= instance )
        reduceP_B_I[name].resize(instance+1);
    if( reduceP_C_I[name].size() <= instance )
        reduceP_C_I[name].resize(instance+1);
    if( reduceP_D_I[name].size() <= instance )
        reduceP_D_I[name].resize(instance+1);
    if( reduceFecundity[name].size() <= instance )
        reduceFecundity[name].resize(instance+1);

    if( elt.getAvailabilityReduction().present() ){
        const scnXml::AvailabilityReduction& elt2 = elt.getAvailabilityReduction().get();
        if( elt2.getInitial() > 1.0 )
            throw util::xml_scenario_error( "availabilityReduction intervention: initial effect must be <= 1" );
        reduceNHHAvailability[name][instance].set (elt2.getInitial(), decay, "availabilityReduction");
    }
    if( elt.getPreprandialKillingEffect().present() ){
        const scnXml::PreprandialKillingEffect& elt2 = elt.getPreprandialKillingEffect().get();
        if( elt2.getInitial() < -1.0 )
            throw util::xml_scenario_error( "reduceAvailability intervention: initial effect must be ≥ -1" );
        reduceP_B_I[name].resize(instance+1);
        reduceP_B_I[name][instance].set (elt2.getInitial(), decay, "reduceP_B_I");
    }
    if( elt.getPostprandialKillingEffect().present() ){
        const scnXml::PostprandialKillingEffect& elt2 = elt.getPostprandialKillingEffect().get();
        if( elt2.getInitial() < -1.0 )
            throw util::xml_scenario_error( "reduceAvailability intervention: initial effect must be ≥ -1" );
        reduceP_C_I[name].resize(instance+1);
        reduceP_C_I[name][instance].set (elt2.getInitial(), decay, "reduceP_C_I");
    }
    if( elt.getRestingKillingEffect().present() ){
        const scnXml::RestingKillingEffect& elt2 = elt.getRestingKillingEffect().get();
        if( elt2.getInitial() < -1.0 )
            throw util::xml_scenario_error( "reduceAvailability intervention: initial effect must be ≥ -1" );
        reduceP_D_I[name].resize(instance+1);
        reduceP_D_I[name][instance].set (elt2.getInitial(), decay, "reduceP_D_I");
    }
    if( elt.getFecundityReduction().present() ){
        const scnXml::FecundityReduction& elt2 = elt.getFecundityReduction().get();
        if( elt2.getInitial() < 0 ||  elt2.getInitial() > 1)
            throw util::xml_scenario_error( "FecundityReduction intervention: initial effect must be be between 0 and 1" );
        reduceFecundity[name][instance].set (elt2.getInitial(), decay, "reduceFecundity");
    }
}
void AnophelesModel::initAddNonHumanHostsInterv(const scnXml::NonHumanHostsVectorSpecies& elt, string name ){
    // Check that the nonHumanHostsType does not exist
    if(addedNhh.count(name) != 0 || initNhh.count(name) != 0)
        throw util::xml_scenario_error( "non human hosts type already exists" );

    NHHParams nhh;
    nhh.mosqRelativeAvailabilityHuman = elt.getMosqRelativeAvailabilityHuman().getValue();
    nhh.mosqProbBiting = elt.getMosqProbBiting().getValue();
    nhh.mosqProbFindingRestSite = elt.getMosqProbFindRestSite().getValue();
    nhh.mosqProbResting = elt.getMosqProbResting().getValue();
    addedNhh[name] = nhh;
}

void AnophelesModel::deployVectorPopInterv (LocalRng& rng, size_t instance){
    assert( instance < emergenceReduction.size() );
    emergenceReduction[instance].deploy( rng, sim::now() );
    // do same as in above function (of EmergenceModel)
    assert( instance < seekingDeathRateIntervs.size() && instance < probDeathOvipositingIntervs.size() );
    seekingDeathRateIntervs[instance].deploy( rng, sim::now() );
    probDeathOvipositingIntervs[instance].deploy( rng, sim::now() );
}
void AnophelesModel::deployVectorTrap(LocalRng& rng, size_t species, size_t instance, double popSize, SimTime lifespan){
    assert(instance < trapParams.size());
    TrapData data;
    data.instance = instance;
    double adultAvail = PerHostAnophParams::get(species).entoAvailability.mean();
    data.initialAvail = popSize * adultAvail * trapParams[instance].relAvail;
    data.availHet = trapParams[instance].availDecay->hetSample(rng);
    data.deployTime = sim::now();
    data.expiry = sim::now() + lifespan;
    baitedTraps.push_back(data);
}
void AnophelesModel::deployNonHumanHostsInterv(LocalRng& rng, size_t species, size_t instance, string name){
    if(initNhh.count(name) == 0)
        throw util::xml_scenario_error("non human hosts type "+name+" not deployed during non human hosts intervention deployment");

    reduceNHHAvailability[name][instance].deploy( rng, sim::now() );
    reduceP_B_I[name][instance].deploy( rng, sim::now() );
    reduceP_C_I[name][instance].deploy( rng, sim::now() );
    reduceP_D_I[name][instance].deploy( rng, sim::now() );
    reduceFecundity[name][instance].deploy( rng, sim::now() );
}

void AnophelesModel::deployAddNonHumanHosts(LocalRng& rng, size_t species, string name, double popSize, SimTime lifespan){
    if(initNhh.count(name) != 0)
        throw util::xml_scenario_error("non human hosts type "+name+" already deployed during non human hosts deployment");

    const NHHParams &nhhParams = addedNhh[name];

    double adultAvail = PerHostAnophParams::get(species).entoAvailability.mean();
    double avail_i = popSize * adultAvail * nhhParams.mosqRelativeAvailabilityHuman;

    NHH nhh;
    nhh.avail_i = avail_i;
    nhh.P_B_I = nhhParams.mosqProbBiting;
    nhh.P_C_I = nhhParams.mosqProbFindingRestSite;
    nhh.P_D_I = nhhParams.mosqProbResting;
    nhh.rel_fecundity = 1.0;
    nhh.expiry = sim::now() + lifespan;
    initNhh[name] = nhh;
}
// Every SimTime::oneTS() days:
void AnophelesModel::advancePeriod (
        double sum_avail, double sigma_df, vector<double>& sigma_dif, double sigma_dff, bool isDynamic)
{
    interventionSurvival = 1.0;
    for( size_t i = 0; i < emergenceReduction.size(); ++i ){
        interventionSurvival *= 1.0 - emergenceReduction[i].current_value( sim::ts0() );
    }
    
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
    for( const util::SimpleDecayingValue& increase : seekingDeathRateIntervs ){
        leaveRate *= 1.0 + increase.current_value( sim::ts0() );
    }
    leaveRate += sum_avail;

    // NON-HUMAN HOSTS INTERVENTIONS
    // Check if some nhh must be removed
    for( auto it = initNhh.begin(); it != initNhh.end();){
        if( sim::ts0() >= it->second.expiry ){
            it = initNhh.erase(it);
            continue;
        }
        it++;
    }

    double modified_nhh_avail = 0.0;
    double modified_nhh_sigma_df = 0.0;
    double modified_nhh_sigma_dff = 0.0;

    map<string,NHH> currentNhh = initNhh;

    for( auto it = reduceNHHAvailability.begin(); it != reduceNHHAvailability.end(); ++it) {
        for( const auto &decay : it->second )
        {
            if(currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].avail_i *= 1.0 - decay.current_value( sim::ts0() );
        }
    }

    for( auto it = reduceP_B_I.begin(); it != reduceP_B_I.end(); ++it) {
        for( const auto &decay : it->second )
        {
            if(currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].P_B_I *= 1.0 - decay.current_value( sim::ts0() );
        }
    }

    for( auto it = reduceP_C_I.begin(); it != reduceP_C_I.end(); ++it) {
        for( const auto &decay : it->second )
        {
            if(currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].P_C_I *= 1.0 - decay.current_value( sim::ts0() );
        }
    }

    for( auto it = reduceP_D_I.begin(); it != reduceP_D_I.end(); ++it) {
        for( const auto &decay : it->second )
        {
            if(currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].P_D_I *= 1.0 - decay.current_value( sim::ts0() );
        }
    }

    for( auto it = reduceFecundity.begin(); it != reduceFecundity.end(); ++it) {
        for( const auto &decay : it->second )
        {
            if(currentNhh.count(it->first) != 0) // Check that the non-human hosts still exist
                currentNhh[it->first].rel_fecundity *= 1.0 - decay.current_value( sim::ts0() );
        }
    }

    for( auto it = currentNhh.begin(); it != currentNhh.end(); ++it){
        modified_nhh_avail += it->second.avail_i;
        const double df = it->second.avail_i * it->second.P_B_I * it->second.P_C_I * it->second.P_D_I;    // term in P_df series
        modified_nhh_sigma_df += df;
        modified_nhh_sigma_dff += df * it->second.rel_fecundity;
    }

    leaveRate += modified_nhh_avail;
    sigma_df += modified_nhh_sigma_df;
    sigma_dff += modified_nhh_sigma_dff;
    // NON-HUMAN HOSTS INTERVENTIONS

    for( auto it = baitedTraps.begin(); it != baitedTraps.end();){
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
    for( const util::SimpleDecayingValue& pDeath : probDeathOvipositingIntervs ){
        alphaE *= 1.0 - pDeath.current_value( sim::ts0() );
    }
    double tsP_df  = sigma_df * alphaE;
    double tsP_dff = sigma_dff * alphaE;

    // from now, sigma_dif becomes P_dif (but we can't simply rename):
    vectors::scale( sigma_dif, alphaE );
    
    // Summed per day:
    partialEIR.assign( WithinHost::Genotypes::N(), 0.0 );
    
    resetTSStats();
    
    // Computing for output only
    double tsP_Amu = (1-tsP_A) * mosqSeekingDeathRate/(mosqSeekingDeathRate + sum_avail + modified_nhh_avail);
    double tsP_A1 = (1-tsP_A) * sum_avail/(mosqSeekingDeathRate + sum_avail + modified_nhh_avail);
    double tsP_Ah = 0.0;
    for( auto it = currentNhh.begin(); it != currentNhh.end(); ++it){
        tsP_Ah += (1-tsP_A) * it->second.avail_i / (mosqSeekingDeathRate + sum_avail + modified_nhh_avail);
    }

    // The code within the for loop needs to run per-day, wheras the main
    // simulation uses one or five day time steps.
    const SimTime nextTS = sim::ts0() + SimTime::oneTS();
    for( SimTime d0 = sim::ts0(); d0 < nextTS; d0 += SimTime::oneDay() ){
        update( d0, tsP_A, tsP_Amu, tsP_A1, tsP_Ah, tsP_df, sigma_dif, tsP_dff, isDynamic, partialEIR, availDivisor);
    }
}

void AnophelesModel::update( SimTime d0, double tsP_A, double tsP_Amu, double tsP_A1, double tsP_Ah, double tsP_df,
        const vector<double> tsP_dif, double tsP_dff,
        bool isDynamic,
        vector<double>& partialEIR, double EIR_factor)
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
    P_Amu[t1] = tsP_Amu;
    P_A1[t1] = tsP_A1;
    P_Ah[t1] = tsP_Ah;
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

            // We use time at end of step (i.e. start + 1) in index:
    SimTime d5Year = mod_nn(d1, SimTime::fromYearsI(5));
    quinquennialS_v[d5Year] = total_S_v;

    const double nOvipositing = P_dff[ttau] * N_v[ttau];       // number ovipositing on this step
    const double newAdults = emergence->update(d0, mosqEmergeRate, nOvipositing) * interventionSurvival;
    util::streamValidate( newAdults );

    // num seeking mosquitos is: new adults + those which didn't find a host
    // yesterday + those who found a host tau days ago and survived cycle:
    N_v[t1] = newAdults + P_A[t0]  * N_v[t0] + nOvipositing;
    
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

void AnophelesModel::uninfectVectors() {
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
double AnophelesModel::getLastVecStat( VecStat vs )const{
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
        case PAmu: return sum1(P_Amu, end, N_v_length);
        case PA1: return sum1(P_A1, end, N_v_length);
        case PAh: return sum1(P_Ah, end, N_v_length);
        default: throw SWITCH_DEFAULT_EXCEPTION;
    }
}
void AnophelesModel::summarize( size_t species )const{
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
