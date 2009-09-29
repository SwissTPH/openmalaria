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
#include "Transmission/Vector/VectorEmergence.h"
#include "Transmission/PerHostTransmission.h"
#include "inputData.h"
#include "human.h"
#include "simulation.h"
#include "util/vectors.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include <fstream>
#include <ctime>


string VectorAnopheles::initialise (const scnXml::Anopheles& anoph, size_t sIndex, const std::list<Human>& population, int populationSize, vector<double>& initialisationEIR) {
  // -----  Set model variables  -----
  const scnXml::Mosq mosq = anoph.getMosq();
  
  mosqRestDuration = mosq.getMosqRestDuration();
  EIPDuration = mosq.getExtrinsicIncubationPeriod();
  
  mosqSeekingDeathRate = mosq.getMosqSeekingDeathRate();
  mosqSeekingDuration = mosq.getMosqSeekingDuration();
  
  humanBase = mosq;
  probMosqSurvivalOvipositing = mosq.getMosqProbOvipositing();
  
  if (1 > mosqRestDuration || mosqRestDuration > EIPDuration) {
    throw xml_scenario_error ("Code expects EIPDuration >= mosqRestDuration >= 1");
  }
  N_v_length = EIPDuration + mosqRestDuration;
  
  const scnXml::Anopheles::NonHumanHostsSequence& otherHosts = anoph.getNonHumanHosts();
  nonHumanHosts.resize (otherHosts.size());
  for (size_t i = 0; i < otherHosts.size(); ++i) {
    nonHumanHosts[i] = otherHosts[i];
  }
  
  
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
  
  
  // -----  if we're driving initialisation from EIR data  -----
  if (anoph.getEir().present()) {
    const scnXml::Eir& eirData = anoph.getEir().get();
    FCEIR.resize (5);
    FCEIR[0] = eirData.getA0();
    FCEIR[1] = eirData.getA1();
    FCEIR[2] = eirData.getB1();
    FCEIR[3] = eirData.getA2();
    FCEIR[4] = eirData.getB2();
    EIRRotateAngle = eirData.getEIRRotateAngle();
    
    vector<double> speciesEIR (Global::intervalsPerYear);
    
    // Calculate forced EIR for pre-intervention phase from FCEIR:
    calcInverseDFTExp(speciesEIR, FCEIR);
    
    if(EIRRotateAngle != 0.0)
      rotateArray(speciesEIR, EIRRotateAngle);
    
    // Add to the TransmissionModel's EIR, used for the initalization phase:
    for (size_t i = 0; i < Global::intervalsPerYear; ++i)
      initialisationEIR[i] += speciesEIR[i];
  }
  
  
  // -----  if we have emerge rate data to use or validate  -----
  if (anoph.getEmergence().present()) {
    const scnXml::Emergence& emergeData = anoph.getEmergence().get();
    
    mosqEmergeRate = vectors::DoubleList2std (emergeData.getEmergenceRate(), daysInYear);
    vectors::scale (mosqEmergeRate, populationSize);
    
    P_dif = vectors::DoubleList2std (emergeData.getKappa(), N_v_length);
    if (!FCEIR.size()) {
      cerr << "Warning: running the vector code during initialisation probably isn't valid." << endl;
      /* Explanation:
       * Normally an entomological infection rate of mosquitoes to humans is
       * forced during initialisation. Without this, some infections are
       * introduced into humans by initially infected mosquitoes (O_v and S_v),
       * but since humans are initially uninfected and have no immunity, no new
       * mosquitoes will be infected initially. There has been some argument
       * that the whole system should in any case reach an equilibrium state by
       * the end of the initialisation phase, but there are no guarantees since
       * the whole system is dynamic, not forced.
       * 
       * Currently this is used for a unittest. */
      initFeedingCycleProbs (sIndex, population, P_dif);
    }
    //else: kappa is validated (from P_dif) and P_* calculated by initMainSimulation
    
    N_v = vectors::DoubleList2std (emergeData.getN_v(), N_v_length);
    vectors::scale (N_v, populationSize);
    O_v = vectors::DoubleList2std (emergeData.getO_v(), N_v_length);
    vectors::scale (O_v, populationSize);
    S_v = vectors::DoubleList2std (emergeData.getS_v(), N_v_length);
    vectors::scale (S_v, populationSize);
  }
  
  return anoph.getMosquito();
}

void VectorAnopheles::destroy () {
}

double VectorAnopheles::calcCycleProbabilities (double& intP_A, double& intP_df, double& intP_dif, size_t sIndex, const std::list<Human>& population, bool withoutInfectiousnessP_dif = false) {
  // rate at which mosquitoes find hosts or die (i.e. leave host-seeking state)
  double leaveSeekingStateRate = mosqSeekingDeathRate;
  for (std::list<Human>::const_iterator h = population.begin(); h != population.end(); ++h)
    leaveSeekingStateRate += h->perHostTransmission.entoAvailability(humanBase, sIndex, h->getAgeInYears());
  for (NonHumanHostsType::const_iterator nnh = nonHumanHosts.begin(); nnh != nonHumanHosts.end(); ++nnh)
    leaveSeekingStateRate += nnh->entoAvailability;
  
  // Probability of a mosquito not finding a host this day:
  intP_A = exp(-leaveSeekingStateRate * mosqSeekingDuration);
  
  // NC's non-autonomous model provides two methods for calculating P_df and
  // P_dif; here we assume that P_E is constant.
  intP_df = 0.0;
  intP_dif = 0.0;
  for (std::list<Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
    const PerHostTransmission& host = h->perHostTransmission;
    double prod = host.entoAvailability(humanBase, sIndex, h->getAgeInYears()) *
      host.probMosqBiting(humanBase, sIndex) *
      host.probMosqResting(humanBase, sIndex);
    intP_df += prod;
    intP_dif += prod * h->probTransmissionToMosquito();
  }
  
  /* For initialisation, we multiply later by the whole population's kappa.
  In this case this function is only called once so performance is less
  important (this allows avoiding code duplication). */
  if (withoutInfectiousnessP_dif)
    intP_dif = intP_df;
  
  for (NonHumanHostsType::const_iterator nnh = nonHumanHosts.begin(); nnh != nonHumanHosts.end(); ++nnh) {
    intP_df += nnh->entoAvailability * nnh->probMosqBitingAndResting();
    // Note: in model, we do the same for intP_dif, except in this case it's
    // multiplied by infectiousness of host to mosquito which is zero.
  }
  
  double P_Ai_base = (1.0 - intP_A) / leaveSeekingStateRate;
  intP_df  *= P_Ai_base * probMosqSurvivalOvipositing;
  intP_dif *= P_Ai_base * probMosqSurvivalOvipositing;
  return P_Ai_base;
}

void VectorAnopheles::initFeedingCycleProbs (size_t sIndex, const std::list<Human>& population, vector<double>& kappaDaily) {
  //Note: kappaDaily may be [a reference to] P_dif, so don't access it after setting P_dif.
  
  double intP_A, intP_df, intP_dif;
  calcCycleProbabilities (intP_A, intP_df, intP_dif, sIndex, population, true);
  
  P_A  .resize (N_v_length);
  P_df .resize (N_v_length);
  P_dif.resize (N_v_length);
  for (int t = 0; t < N_v_length; ++t) {
    P_A[t]	= intP_A;
    P_df[t]	= intP_df;
    P_dif[t]	= intP_dif * kappaDaily[t];	// FIXME: don't use kappa for initialisation (only use for non-vector model) - initialise by calling advancePeriod or similar per timestep?
  }
}

void VectorAnopheles::initMainSimulation (size_t sIndex, const std::list<Human>& population, int populationSize, vector<double>& kappa) {
  // If EIR data was provided, validate EIR or do emergence rate calculations
  // and switch to calculating a dynamic EIR.
  if (FCEIR.size()) {
    //BEGIN Validate kappa, initialise P_A, P_df, P_dif
    if (Simulation::simulationTime*Global::interval < N_v_length)
      throw xml_scenario_error ("Initialization phase too short");
    vector<double> kappaDaily (N_v_length, 0.0);
    
    // For day over N_v_length days prior to the next timestep's day.
    int endDay = (Simulation::simulationTime+1) * Global::interval;
    for (int day = endDay - N_v_length; day < endDay; ++day) {
      // Should correspond to index of kappa updated by updateKappa:
      kappaDaily[day % N_v_length] =
	kappa[(day / Global::interval - 1) % Global::intervalsPerYear];
    }
    
    bool valid = vectors::approxEqual (kappaDaily, P_dif);
    
    initFeedingCycleProbs (sIndex, population, kappaDaily);
    //END Validate kappa, initialise P_A, P_df, P_dif
    
    /* This (timely) check of emergence rates and N/O/S_v is skipped if
     * a -v or --noErcValidation option is present
     * AND kappaDaily was valid. */
    if (!valid || !(Global::clOptions & CLO::NO_ERC_VALIDATION)) {
      time_t seconds = time (NULL);
      cout << "Starting emergence rate validation" << endl;
      //BEGIN Validate or calculate mosqEmergeRate
      double sumAvailability = 0.0;
      for (std::list<Human>::const_iterator h = population.begin(); h != population.end(); ++h)
	sumAvailability += h->perHostTransmission.entoAvailability(humanBase, sIndex, h->getAgeInYears());
      // Use if you want to recalculate the constant in HostMosquitoInteraction::initialise:
      //cout << "Average human availability: " << setprecision(12) << sumAvailability/populationSize << endl;
      //NOTE: could check averageAvailability is reasonably close to humanBase.entoAvailability
      // (though it's not clear what to do and how much of a problem it is if this isn't the case).
      
      // A class encapsulating VectorInternal code. Destructor frees memory at end of this function.
      VectorEmergence emerge (mosqRestDuration, EIPDuration,
			      populationSize,
			      mosqSeekingDeathRate, mosqSeekingDuration,
			      sumAvailability,
			      humanBase, probMosqSurvivalOvipositing);
      
      
      /* We can either take the EIR from the initialisation stage and expand to
       * length daysInYear without smoothing:
       *	EIRInit = convertLengthToFullYear(speciesEIR);
       * or recalculate from fourier coefficients to get a smooth array: */
      vector<double> EIRInit(daysInYear, 0.0);
      calcInverseDFTExp(EIRInit, FCEIR);
      if(EIRRotateAngle != 0.0)
	rotateArray(EIRInit, EIRRotateAngle);
      
      
      if (mosqEmergeRate.size() == 0) {	// not set; we'll need an estimate
	mosqEmergeRate.resize (daysInYear);
	// The root finding algorithm needs some estimate to start from. It's accuracy isn't that important since it should converge in a single step anyway.
	double temp = populationSize*sumAvailability;
	for (int i = 0; i < daysInYear; i++) {
	  mosqEmergeRate[i] = EIRInit[i]*temp;
	}
      }
      
      //TODO validate as separate function, do all validation first, then all calculations if anything isn't valid
      
      if (!emerge.CalcInitMosqEmergeRate(convertLengthToFullYear (kappa),
					 EIRInit,
					 nonHumanHosts,
					 mosqEmergeRate)) {
	valid = false;	// test failed; emergence rates recalculated
      }
      //END Validate or calculate mosqEmergeRate
      
      //BEGIN Get and validate N_v, O_v and S_v
      if (Simulation::simulationTime*Global::interval < N_v_length || N_v_length > daysInYear)
	throw xml_scenario_error ("Initialization phase or daysInYear too short");
      vector<double> Nv(N_v_length), Ov(N_v_length), Sv(N_v_length);
      
      // Retrieve the periodic orbits for Nv, Ov, and Sv.
      gsl_vector** x_p;
      size_t mt = emerge.getN_vO_vS_v (x_p);
      for (int day = endDay - N_v_length; day < endDay; ++day) {
	size_t t = day % N_v_length;
	size_t i = (daysInYear + day - endDay) % daysInYear;
	
	Nv[t] = gsl_vector_get(x_p[i], 0);
	Ov[t] = gsl_vector_get(x_p[i], mt);
	Sv[t] = gsl_vector_get(x_p[i], 2*mt);
      }
      
      valid = valid
	&& vectors::approxEqual (N_v, Nv)
	&& vectors::approxEqual (O_v, Ov)
	&& vectors::approxEqual (S_v, Sv);
      
      if (N_v.size()) {
	if (valid)
	  cerr << "Emergence rate parameters in scenario document were accurate." << endl;
	else
	  cerr << "Warning: emergence rate parameters in scenario document were not accurate." << endl;
      }
      
      // Set whether or not valid; no harm if they already have good values
      N_v = Nv;
      O_v = Ov;
      S_v = Sv;
      //END Get and validate N_v, O_v and S_v
      
      //BEGIN Write out new values
      if (!valid) {
	cerr << "Parameters associated with emergence rate have been generated will be saved." << endl;;
	
	scnXml::DoubleList sxEmergeRate, sxKappa, sxNv, sxOv, sxSv;
	typedef scnXml::DoubleList::ItemSequence DLIS;
	
	vector<double> tp (mosqEmergeRate);
	vectors::scale (tp, 1.0 / populationSize);
	sxEmergeRate.setItem (DLIS (tp.begin(), tp.end()));
	
	sxKappa.setItem (DLIS (kappaDaily.begin(), kappaDaily.end()));
	
	tp = Nv;
	vectors::scale (tp, 1.0 / populationSize);
	sxNv.setItem (DLIS (tp.begin(), tp.end()));
	
	tp = Ov;
	vectors::scale (tp, 1.0 / populationSize);
	sxOv.setItem (DLIS (tp.begin(), tp.end()));
	
	tp = Sv;
	vectors::scale (tp, 1.0 / populationSize);
	sxSv.setItem (DLIS (tp.begin(), tp.end()));
	
	scnXml::Emergence sxEmerge (sxEmergeRate, sxKappa, sxNv, sxOv, sxSv);
	
	// Assumptions: EntoData.Vector element exists, anophleses[sIndex] exists
	scnXml::Anopheles& sxAnoph = getMutableScenario()
	  .getEntoData()
	  .getVector().get()
	  .getAnopheles()[sIndex];
	sxAnoph.setEmergence (sxEmerge);
	documentChanged = true;
      }
      //END Write out new values
      cout << "Finished validating/calculation emergence rates ("<< (time (NULL) - seconds) <<" seconds)" << endl;
    }
    
    FCEIR.resize(0, 0.0);	// indicate EIR-driven mode is no longer being used (although this is ignored)
  }
  // else: we're already in dynamic EIR calculation mode
}


// Every Global::interval days:
void VectorAnopheles::advancePeriod (const std::list<Human>& population, int simulationTime, size_t sIndex) {
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
  Let P_Ai_base[t] = (1 - P_A[t]) / sum_{h in hosts} α_h[t].
  
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
  double intP_A, intP_df, intP_dif;
  double P_Ai_base = calcCycleProbabilities (intP_A, intP_df, intP_dif, sIndex, population);
  //END P_A, P_Ai, P_df, P_dif
  
  // Summed per day:
  partialEIR = 0.0;
  
  
  // The code within the for loop needs to run per-day, wheras the main
  // simulation uses Global::interval day (currently 5 day) time steps.
  int endDay = (simulationTime+1) * Global::interval;
  for (int day = simulationTime * Global::interval; day < endDay; ++day) {
    // Indecies for today, yesterday and mosqRestDuration days back:
    // Warning: with x<0, x%y can be negative (depending on compiler); avoid x<0.
    // Use dMod rather than day when subtracting something then using modulus
    size_t dMod = day + N_v_length;
    assert (dMod >= (size_t)N_v_length);
    size_t t    = dMod % N_v_length;
    size_t t1   = (dMod - 1) % N_v_length;
    size_t ttau = (dMod - mosqRestDuration) % N_v_length;
    
    
    // These only need to be calculated once per timestep, but should be
    // present in each of the previous N_v_length - 1 positions of arrays.
    P_A[t] = intP_A;
    P_df[t] = intP_df;
    P_dif[t] = intP_dif;
    
    
    N_v[t] = mosqEmergeRate[day%daysInYear] * larvicidingIneffectiveness
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
    
    
    partialEIR += S_v[t] * P_Ai_base;
  }
}


void VectorAnopheles::intervLarviciding (const scnXml::LarvicidingAnopheles& elt) {
  larvicidingIneffectiveness = 1 - elt.getEffectiveness();
  larvicidingEndStep = Simulation::simulationTime + (elt.getDuration() / Global::interval);
}


vector<double> VectorAnopheles::convertLengthToFullYear (vector<double>& ShortArray) {
  vector<double> FullArray (daysInYear);
  for (size_t i=0; i < Global::intervalsPerYear; i++) {
    for (int j=0; j < Global::interval; j++) {
      FullArray[i*Global::interval + j] = ShortArray[i];
    }
  }
  return FullArray;
}


void VectorAnopheles::calcInverseDFTExp(vector<double>& tArray, vector<double>& FC) {
  if (FC.size() % 2 == 0)
    throw xml_scenario_error("The number of Fourier coefficents should be odd.");
  
  // Frequency
  double w = 2*M_PI / double(tArray.size());
  
  // Number of Fourier Modes.
  int Fn = (FC.size()-1)/2;
  
  // Calculate inverse discrete Fourier transform
  for (size_t t=0; t<tArray.size(); t++){
    double temp = FC[0];
    double wt = w*(t+1);
    for(int n=1;n<=Fn;n++){
      temp = temp + FC[2*n-1]*cos(n*wt) + FC[2*n]*sin(n*wt);
    }
    tArray[t] = exp(temp);
  }
}

void VectorAnopheles::rotateArray(vector<double>& rArray, double rAngle) {
  vector<double> tempArray (rArray.size());
  size_t rotIndex = size_t ((rAngle*rArray.size())/(2.0*M_PI));

  for (size_t i=0; i < rArray.size(); i++) {
    tempArray[(i+rotIndex) % rArray.size()] = rArray[i];
  }
  rArray = tempArray;
}
