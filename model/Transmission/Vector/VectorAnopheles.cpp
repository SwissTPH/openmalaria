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
#include "Host/Human.h"

#include "inputData.h"
#include "Simulation.h"
#include "util/vectors.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include <fstream>
#include <ctime>

namespace OM { namespace Transmission {
    using namespace OM::util;
    
string VectorAnopheles::initialise (const scnXml::Anopheles& anoph, size_t sIndex, vector<double>& initialisationEIR) {
  // -----  Set model variables  -----
  const scnXml::Mosq mosq = anoph.getMosq();
  
  mosqRestDuration = mosq.getMosqRestDuration();
  EIPDuration = mosq.getExtrinsicIncubationPeriod();
  
  mosqSeekingDeathRate = mosq.getMosqSeekingDeathRate();
  mosqSeekingDuration = mosq.getMosqSeekingDuration();
  
  humanBase = mosq;
  probMosqSurvivalOvipositing = mosq.getMosqProbOvipositing();
  
  if (1 > mosqRestDuration || mosqRestDuration > EIPDuration) {
    throw util::xml_scenario_error ("Code expects EIPDuration >= mosqRestDuration >= 1");
  }
  N_v_length = EIPDuration + mosqRestDuration;
  
  const scnXml::Anopheles::NonHumanHostsSequence& otherHosts = anoph.getNonHumanHosts();
  nonHumanHosts.resize (otherHosts.size());
  for (size_t i = 0; i < otherHosts.size(); ++i) {
    nonHumanHosts[i] = otherHosts[i];
  }
  
  
  // -----  EIR  -----
  const scnXml::Eir& eirData = anoph.getEir();
  
  FSCoeffic.resize (5);
  FSCoeffic[0] = eirData.getA0();
  FSCoeffic[1] = eirData.getA1();
  FSCoeffic[2] = eirData.getB1();
  FSCoeffic[3] = eirData.getA2();
  FSCoeffic[4] = eirData.getB2();
  EIRRotateAngle = eirData.getEIRRotateAngle();
  
  // Calculate forced EIR for pre-intervention phase from FSCoeffic:
  vector<double> speciesEIR (Global::intervalsPerYear);
  calcFourierEIR (speciesEIR, FSCoeffic, EIRRotateAngle);
  vectors::scale (speciesEIR, Global::interval);	// input EIR is per-capita per-day, so scale to per-interval
  
  // Add to the TransmissionModel's EIR, used for the initalization phase:
  for (size_t i = 0; i < Global::intervalsPerYear; ++i)
    initialisationEIR[i] += speciesEIR[i];
  
  // Set other data used for mosqEmergeRate calculation:
  FSRotateAngle = EIRRotateAngle - (EIPDuration+10)/365.*2.*M_PI;	// usually around 20 days; no real analysis for effect of changing EIPDuration or mosqRestDuration
  initNvFromSv = 1.0 / anoph.getPropInfectious();
  initNv0FromSv = initNvFromSv * anoph.getPropInfected();	// temporarily use of initNv0FromSv
  
  
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
  
  annualS_v.assign (Global::DAYS_IN_YEAR, 0.0);
  forcedS_v.resize (Global::DAYS_IN_YEAR);
  mosqEmergeRate.resize (Global::DAYS_IN_YEAR);	// Only needs to be done here if loading from checkpoint
  
  
  return anoph.getMosquito();
}
void VectorAnopheles::setupNv0 (size_t sIndex, const std::list<Host::Human>& population, int populationSize) {
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
    double prod = host.entoAvailabilityFull (humanBase, sIndex, h->getAgeInYears(), transmissionModel->ageCorrectionFactor);
    leaveSeekingStateRate += prod;
    prod *= host.probMosqBiting(humanBase, sIndex);
    sumPFindBite += prod;
    intP_df += prod * host.probMosqResting(humanBase, sIndex);
  }
  
  for (NonHumanHostsType::const_iterator nnh = nonHumanHosts.begin(); nnh != nonHumanHosts.end(); ++nnh) {
    leaveSeekingStateRate += nnh->entoAvailability;
    intP_df += nnh->entoAvailability * nnh->probMosqBitingAndResting();
    // Note: in model, we do the same for intP_dif, except in this case it's
    // multiplied by infectiousness of host to mosquito which is zero.
  }
  
  // Probability of a mosquito not finding a host this day:
  double intP_A = exp(-leaveSeekingStateRate * mosqSeekingDuration);
  
  double P_Ai_base = (1.0 - intP_A) / leaveSeekingStateRate;
  sumPFindBite *= P_Ai_base;
  intP_df  *= P_Ai_base * probMosqSurvivalOvipositing;
  //END P_A, P_Ai, P_df, P_dif
  
  double initOvFromSv = initNv0FromSv;	// temporarily use of initNv0FromSv
  
  FSCoeffic[0] += log (populationSize / sumPFindBite);	// same as multiplying resultant eir since calcFourierEIR takes exp(...)
  calcFourierEIR (forcedS_v, FSCoeffic, FSRotateAngle);
  
  mosqEmergeRate = forcedS_v;
  initNv0FromSv = initNvFromSv * (1.0 - intP_A - intP_df);
  vectors::scale (mosqEmergeRate, initNv0FromSv);
  
  for (int t = 0; t < N_v_length; ++t) {
    P_A[t] = intP_A;
    P_df[t] = intP_df;
    P_dif[t] = 0.0;	// humans start off with no infectiousness.. so just wait
    //t in: sTime*Global::interval..((sTime+1)*Global::interval-1)
    int sTime = t / Global::interval;
    S_v[t] = forcedS_v[sTime % Global::intervalsPerYear];
    N_v[t] = S_v[t] * initNvFromSv;
    O_v[t] = S_v[t] * initOvFromSv;
  }
  
  sumAnnualForcedS_v = vectors::sum (forcedS_v);
  
  // All set up to drive simulation from forcedS_v
}

void VectorAnopheles::destroy () {
}

bool VectorAnopheles::vectorInitIterate () {
  // Try to match S_v against its predicted value. Don't try with N_v or O_v
  // because the predictions will change - would be chasing a moving target!
  // EIR comes directly from S_v, so should fit after we're done.
  
  double factor = sumAnnualForcedS_v / vectors::sum(annualS_v);
//   cout << "Pre-calced Sv, dynamic Sv:\t"<<sumAnnualForcedS_v<<'\t'<<vectors::sum(annualS_v)<<endl;
  if (!(factor > 1e-6 && factor < 1e6))	// unlikely, but might as well check incase either operand was zero
    throw runtime_error ("factor out of bounds");
  
  //NOTE: how accurate do we want fits to be? Results seem to be stocastically something like 1-5% different.
  const double LIMIT = 0.05;
  if (fabs(factor - 1.0) > LIMIT) {
    cout << "Vector iteration: adjusting with factor "<<factor<<endl;
    // Adjusting mosqEmergeRate is the important bit. The rest should just bring things to a stable state quicker.
    initNv0FromSv *= factor;
    initNvFromSv *= factor;	//(not currently used)
    vectors::scale (mosqEmergeRate, factor);
    vectors::scale (N_v, factor);
    // What factor exactly these should be scaled by isn't obvious; in any case they should reach stable values quickly.
    //vectors::scale (O_v, factor);
    //vectors::scale (S_v, factor);
    return true;	// we scaled mosqEmergeRate; iterate again
  } else {
    // Once the amplitude is approximately correct, we try to find the offset
    double rAngle = Nv0DelayFitting::fit<double> (EIRRotateAngle, FSCoeffic, annualS_v);
    cout << "Vector iteration: rotating with angle (in radians): " << rAngle << endl;
    FSRotateAngle -= rAngle;	// annualS_v was already rotated by old value of FSRotateAngle, so increment
    calcFourierEIR (forcedS_v, FSCoeffic, FSRotateAngle);
    // We use the stored initXxFromYy calculated from the ideal population age-structure (at init).
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);
    
    sumAnnualForcedS_v = vectors::sum (forcedS_v) * Global::interval;
    
    return (rAngle > LIMIT * 2*M_PI / Global::intervalsPerYear);	// iterate again if result wasn't close
  }
}


// Every Global::interval days:
void VectorAnopheles::advancePeriod (const std::list<Host::Human>& population, int simulationTime, size_t sIndex, bool isDynamic) {
  if (simulationTime >= larvicidingEndStep) {
    larvicidingEndStep = std::numeric_limits<int>::max();
    larvicidingIneffectiveness = 1.0;
  }
  
  
  /* Largely equations correspond to Nakul Chitnis's model in
    "A mathematic model for the dynamics of malaria in
    mosquitoes feeding on a heterogeneous host population" [MMDM]
  section 2, 3.5-3.6, plus extensions to a non-autonomous case from
    "Nonautonomous Difference Equations for Malaria Dynamics
                 in a Mosquito Population" [NDEMD]
  
  We calculate EIR over a 5-day Global::interval as:
    sum_{for t over days} σ_i[t] * s_v[t]
    = sum_... (N_v[t] * P_Ai[t] * P_B_i[t])/(T*N_i[t]) * S_v[t]/N_v[t]
    = sum_... P_Ai[t] * P_B_i[t] * S_v[t]
  (since T == 1 and N_i[t] == 1 for all t).
  
    P_Ai[t] = (1 - P_A[t]) α_i[t] / sum_{h in hosts} α_h[t]
  (letting N_h[t] == 1 for all h,t). The only part of this varying per-host is
    α_i[t] = host.entoAvailability (index, h->getAgeInYears())
    Let P_Ai_base[t] = (1 - P_A[t]) / (sum_{h in hosts} α_h[t] + μ_vA).
  
  Note that although the model allows α_i and P_B_i to vary per-day, they only
  vary per Global::interval of the main simulation. Hence:
    EIR = (sum_{t=...} S_v[t] * P_Ai_base[t]) * α_i * P_B_i
  
  Since S_v[t] * P_Ai_base[t] does not vary per individual, we calculate this
  per Global::interval of the main simulation as partialEIR:
    partialEIR = (sum_{t=...} S_v[t] * P_Ai_base[t])
  
  Hence calculateEIR() only needs to do the following:
    EIR = partialEIR * α_i * P_B_i
  */
  
  
  //BEGIN P_A, P_Ai, P_df, P_dif
  // rate at which mosquitoes find hosts or die (i.e. leave host-seeking state)
  double leaveSeekingStateRate = mosqSeekingDeathRate;
  
  // NC's non-autonomous model provides two methods for calculating P_df and
  // P_dif; here we assume that P_E is constant.
  double intP_df = 0.0;
  double intP_dif = 0.0;
  for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
    const PerHostTransmission& host = h->perHostTransmission;
    double prod = host.entoAvailabilityFull (humanBase, sIndex, h->getAgeInYears(), transmissionModel->ageCorrectionFactor);
    leaveSeekingStateRate += prod;
    prod *= host.probMosqBiting(humanBase, sIndex)
	  * host.probMosqResting(humanBase, sIndex);
    intP_df += prod;
    intP_dif += prod * h->probTransmissionToMosquito();
  }
  
  for (NonHumanHostsType::const_iterator nnh = nonHumanHosts.begin(); nnh != nonHumanHosts.end(); ++nnh) {
    leaveSeekingStateRate += nnh->entoAvailability;
    intP_df += nnh->entoAvailability * nnh->probMosqBitingAndResting();
    // Note: in model, we do the same for intP_dif, except in this case it's
    // multiplied by infectiousness of host to mosquito which is zero.
  }
  
  // Probability of a mosquito not finding a host this day:
  double intP_A = exp(-leaveSeekingStateRate * mosqSeekingDuration);
  
  double P_Ai_base = (1.0 - intP_A) / leaveSeekingStateRate;
  intP_df  *= P_Ai_base * probMosqSurvivalOvipositing;
  intP_dif *= P_Ai_base * probMosqSurvivalOvipositing;
  //END P_A, P_Ai, P_df, P_dif
  
  //cout << "t"<<simulationTime<<":\t"<<'\t'<<intP_A<<'\t'<<P_Ai_base <<'\t'<< intP_df <<'\t'<<intP_dif<<endl;
  
  // Summed per day:
  partialEIR = 0.0;
  
  timestep_N_v0 = 0.0, timestep_N_v = 0.0, timestep_O_v = 0.0, timestep_S_v = 0.0;
  
  // The code within the for loop needs to run per-day, wheras the main
  // simulation uses Global::interval day (currently 5 day) time steps.
  int firstDay = simulationTime * Global::interval;
  for (size_t i = 0; i < (size_t)Global::interval; ++i) {
    // Warning: with x<0, x%y can be negative (depending on compiler); avoid x<0.
    // We add N_v_length so that ((dMod - x) >= 0) for (x <= N_v_length).
    size_t dMod = i + firstDay + N_v_length;
    assert (dMod >= (size_t)N_v_length);
    // Indecies for today, yesterday and mosqRestDuration days back:
    size_t t    = dMod % N_v_length;
    size_t t1   = (dMod - 1) % N_v_length;
    size_t ttau = (dMod - mosqRestDuration) % N_v_length;
    // Day of year:
    size_t dYear = (firstDay + i) % Global::DAYS_IN_YEAR;
    
    
    // These only need to be calculated once per timestep, but should be
    // present in each of the previous N_v_length - 1 positions of arrays.
    P_A[t] = intP_A;
    P_df[t] = intP_df;
    P_dif[t] = intP_dif;
    
    
    N_v[t] = mosqEmergeRate[dYear] * larvicidingIneffectiveness
        + P_A[t1]  * N_v[t1]
        + P_df[ttau] * N_v[ttau];
    O_v[t] = P_dif[ttau] * (N_v[ttau] - O_v[ttau])
        + P_A[t1]  * O_v[t1]
        + P_df[ttau] * O_v[ttau];
    
    
    //BEGIN S_v
    // Set up array with n in 1..θ_s−1 for f_τ(dMod-n) (NDEMD eq. 1.7)
    size_t fProdEnd = 2*mosqRestDuration;
    for (size_t n = mosqRestDuration+1; n <= fProdEnd; ++n) {
      ftauArray[n] =
          ftauArray[n-1] * P_A[(dMod-n)%N_v_length];
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
      size_t tsl = (ts - l) % N_v_length;	// index dMod - theta_s - l
      sum += P_dif[tsl] * P_df[ttau] * (N_v[tsl] - O_v[tsl]) * ftauArray[EIPDuration+l-mosqRestDuration];
    }
    
    
    // Set up array with n in 1..θ_s−τ for f(dMod-n) (NDEMD eq. 1.6)
    for (int n = 1; n <= mosqRestDuration; ++n) {
      fArray[n] =
          fArray[n-1] * P_A[(dMod-n)%N_v_length];
    }
    fArray[mosqRestDuration] += P_df[ttau];
    
    fProdEnd = EIPDuration-mosqRestDuration;
    for (size_t n = mosqRestDuration+1; n <= fProdEnd; ++n) {
      size_t tn = (dMod-n)%N_v_length;
      fArray[n] =
          P_df[tn] * fArray[n - mosqRestDuration]
          + P_A[tn] * fArray[n-1];
    }
    
    
    ts = ts % N_v_length;	// index dMod - theta_s
    S_v[t] = P_dif[ts] * fArray[EIPDuration-mosqRestDuration] * (N_v[ts] - O_v[ts])
        + sum
        + P_A[t1]*S_v[t1]
        + P_df[ttau]*S_v[ttau];
    //END S_v
    
    annualS_v[dYear] = S_v[t];
    
    partialEIR += S_v[t] * P_Ai_base;
    
    timestep_N_v += N_v[t];
    timestep_N_v0 += mosqEmergeRate[dYear];
    timestep_O_v += O_v[t];
    timestep_S_v += S_v[t];
  }
  
#ifdef OMV_CSV_REPORTING
  (*csvReporting) << timestep_N_v0/Global::interval << '\t' << timestep_N_v/Global::interval << '\t' << timestep_O_v/Global::interval << '\t' << timestep_S_v/Global::interval << '\t';
#endif
}


void VectorAnopheles::intervLarviciding (const scnXml::LarvicidingAnopheles& elt) {
  cerr << "This larviciding implementation isn't valid (according to NC)." << endl;
  larvicidingIneffectiveness = 1 - elt.getEffectiveness();
  larvicidingEndStep = Global::simulationTime + (elt.getDuration() / Global::interval);
}

void VectorAnopheles::summarize (const string speciesName, Survey& survey) {
    survey.set_Vector_Nv0 (speciesName, timestep_N_v0/Global::interval);
    survey.set_Vector_Nv (speciesName, timestep_N_v/Global::interval);
    survey.set_Vector_Ov (speciesName, timestep_O_v/Global::interval);
    survey.set_Vector_Sv (speciesName, timestep_S_v/Global::interval);
}


void VectorAnopheles::calcFourierEIR (vector<double>& tArray, vector<double>& FC, double rAngle) {
  if (FC.size() % 2 == 0)
      throw util::xml_scenario_error("The number of Fourier coefficents should be odd.");
  
  // Frequency
  double w = 2*M_PI / double(tArray.size());
  
  // Number of Fourier Modes.
  int Fn = (FC.size()-1)/2;
  
  // Calculate inverse discrete Fourier transform
  for (size_t t=0; t<tArray.size(); t++){
    double temp = FC[0];
    double wt = w*t - rAngle;
    for(int n=1;n<=Fn;n++){
      temp = temp + FC[2*n-1]*cos(n*wt) + FC[2*n]*sin(n*wt);
    }
    tArray[t] = exp(temp);
  }
}

} }