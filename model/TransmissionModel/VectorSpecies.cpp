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

#include "TransmissionModel.h"
#include "TransmissionModel/VectorSpecies.h"
#include "inputData.h"
#include "human.h"
#include "TransmissionModel/VectorInternal.h"
#include "simulation.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include <fstream>


void VectorTransmissionSpecies::initialise (const scnXml::Anopheles& anoph, vector<double>& initialisationEIR) {
  scnXml::Mosq mosq = anoph.getMosq();
  
  mosqRestDuration = mosq.getMosqRestDuration();
  EIPDuration = mosq.getExtrinsicIncubationPeriod();
  
  mosqSeekingDeathRate = mosq.getMosqSeekingDeathRate();
  mosqSeekingDuration = mosq.getMosqSeekingDuration();
  
  entoAvailability = mosq.getMosqEntoAvailability();
  probMosqBiting = mosq.getMosqProbBiting();
  probMosqFindRestSite = mosq.getMosqProbFindRestSite();
  probMosqSurvivalResting = mosq.getMosqProbResting();
  probMosqSurvivalOvipositing = mosq.getMosqProbOvipositing();
  
  emergenceRateFilename = mosq.getEmergenceRateFilename();
  
  //TODO: initialize these and perhaps other params properly:
  //K_vi
  //P_vi
  
  if (1 > mosqRestDuration || mosqRestDuration > EIPDuration) {
    throw xml_scenario_error ("Code expects EIPDuration >= mosqRestDuration >= 1");
  }
  
  N_v_length = EIPDuration + mosqRestDuration;
  
  size_t l = N_v_length*sizeof(double);
  P_A = new double[N_v_length * 6];	// allocate memory for all arrays
  P_df	= P_A + l;
  P_dif	= P_df + l;
  N_v	= P_dif + l;
  O_v	= N_v + l;
  S_v	= O_v + l;
  
  // Set up fArray and ftauArray. Each step, all elements not set here are
  // calculated, even if they aren't directly used in the end;
  // however all calculated values are used in calculating the next value.
  fArray.resize(EIPDuration-mosqRestDuration+1);
  fArray[0] = 1.0;
  ftauArray.resize(EIPDuration);
  for (int i = 0; i < mosqRestDuration; ++i)
    ftauArray[i] = 0.0;
  ftauArray[mosqRestDuration] = 1.0;
  
  
  scnXml::Eir eirData = anoph.getEir();
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

void VectorTransmissionSpecies::destroy () {
  delete[] P_A;
}

void VectorTransmissionSpecies::initMainSimulation (size_t sIndex, const std::list<Human>& population, int populationSize, vector<double>& kappa) {
  calMosqEmergeRate (populationSize, kappa);	// initialises N_v, O_v, S_v
  
  
  //BEGIN P_A, P_Ai, P_df, P_dif
  // Per Global::interval (hosts don't update per day):
  double leaveHostRate = mosqSeekingDeathRate;
  for (std::list<Human>::const_iterator h = population.begin(); h != population.end(); ++h)
    leaveHostRate += h->perHostTransmission.entoAvailability(sIndex);
  
  // Probability of a mosquito not finding a host this day:
  double intP_A = exp(-leaveHostRate * mosqSeekingDuration);
  
  double P_Ai_base = (1.0 - intP_A) / leaveHostRate;
  
  // NC's non-autonomous model provides two methods for calculating P_df and
  // P_dif; here we assume that P_E is constant.
  double intP_df = 0.0;
  for (std::list<Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
    const PerHostTransmission& host = h->perHostTransmission;
    double prod = host.entoAvailability(sIndex) * host.probMosqBiting(sIndex)
    * host.probMosqFindRestSite(sIndex) * host.probMosqSurvivalResting(sIndex);
    intP_df += prod;
  }
  intP_df  *= P_Ai_base * probMosqSurvivalOvipositing;
  //END P_A, P_Ai, P_df, P_dif
  
  
  if (Simulation::simulationTime*Global::interval < N_v_length)
    throw xml_scenario_error ("Initialization phase too short");
  // For day over N_v_length days prior to the next timestep's day.
  int endDay = (Simulation::simulationTime+1) * Global::interval;
  for (int day = endDay - N_v_length; day < endDay; ++day) {
    size_t t       = day % N_v_length;
    size_t simStep = day / Global::interval;
    
    P_A[t]	= intP_A;
    P_df[t]	= intP_df;
    // Should correspond to index of kappa updated by updateKappa:
    P_dif[t]	= intP_df * kappa[(simStep-1) % Global::intervalsPerYear];
  }
}


// Every Global::interval days:
void VectorTransmissionSpecies::advancePeriod (const std::list<Human>& population, int simulationTime, size_t sIndex) {
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
    α_i[t] = host.entoAvailability ()
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
  // Per Global::interval (hosts don't update per day):
  double leaveHostRate = mosqSeekingDeathRate;
  for (std::list<Human>::const_iterator h = population.begin(); h != population.end(); ++h)
    leaveHostRate += h->perHostTransmission.entoAvailability(sIndex);
  
  // Probability of a mosquito not finding a host this day:
  double intP_A = exp(-leaveHostRate * mosqSeekingDuration);
  
  double P_Ai_base = (1.0 - intP_A) / leaveHostRate;
  
  // NC's non-autonomous model provides two methods for calculating P_df and
  // P_dif; here we assume that P_E is constant.
  double intP_df = 0.0;
  double intP_dif = 0.0;
  for (std::list<Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
    const PerHostTransmission& host = h->perHostTransmission;
    double prod = host.entoAvailability(sIndex) * host.probMosqBiting(sIndex)
    * host.probMosqFindRestSite(sIndex) * host.probMosqSurvivalResting(sIndex);
    intP_df += prod;
    intP_dif += prod * h->withinHostModel->getProbTransmissionToMosquito();
  }
  intP_df  *= P_Ai_base * probMosqSurvivalOvipositing;
  intP_dif *= P_Ai_base * probMosqSurvivalOvipositing;
  //END P_A, P_Ai, P_df, P_dif
  
  // Summed per day:
  partialEIR = 0.0;
  
  
  // The code within the for loop needs to run per-day, wheras the main
  // simulation uses Global::interval day (currently 5 day) time steps.
  int endDay = (simulationTime+1) * Global::interval;
  for (int day = simulationTime * Global::interval; day < endDay; ++day) {
    // Indecies for today, yesterday and mosqRestDuration days back:
    // Warning: with x<0, x%y can be negative (depending on compiler); avoid x<0.
    size_t t    = day % N_v_length;
    size_t t1   = (day - 1) % N_v_length;
    size_t ttau = (day - mosqRestDuration) % N_v_length;
    
    
    // These only need to be calculated once per timestep, but should be
    // present in each of the previous N_v_length - 1 positions of arrays.
    P_A[t] = intP_A;
    P_df[t] = intP_df;
    P_dif[t] = intP_dif;
    
    
    N_v[t] = mosqEmergeRate[day%daysInYear]
        + P_A[t1]  * N_v[t1]
        + P_df[ttau] * N_v[ttau];
    O_v[t] = P_dif[ttau] * (N_v[ttau] - O_v[ttau])
        + P_A[t1]  * O_v[t1]
        + P_df[ttau] * O_v[ttau];
    
    
    //BEGIN S_v
    // Set up array with n in 1..θ_s−1 for f_τ(day-n) (NDEMD eq. 1.7)
    size_t fProdEnd = 2*mosqRestDuration;
    for (size_t n = mosqRestDuration+1; n <= fProdEnd; ++n) {
      ftauArray[n] =
          ftauArray[n-1] * P_A[(day-n)%N_v_length];
    }
    ftauArray[fProdEnd] += P_df[(day-fProdEnd)%N_v_length];
    
    for (int n = fProdEnd+1; n < EIPDuration; ++n) {
      size_t tn = (day-n)%N_v_length;
      ftauArray[n] =
          P_df[tn] * ftauArray[n - mosqRestDuration]
          + P_A[tn] * ftauArray[n-1];
    }
    
    double sum = 0.0;
    size_t ts = day - EIPDuration;
    for (int l = 1; l < mosqRestDuration; ++l) {
      size_t tsl = (ts - l) % N_v_length;	// index day - theta_s - l
      sum += P_dif[tsl] * P_df[ttau] * (N_v[tsl] - O_v[tsl]) * ftauArray[EIPDuration+l-mosqRestDuration];
    }
    
    
    // Set up array with n in 1..θ_s−τ for f(day-n) (NDEMD eq. 1.6)
    for (int n = 1; n <= mosqRestDuration; ++n) {
      fArray[n] =
          fArray[n-1] * P_A[(day-n)%N_v_length];
    }
    fArray[mosqRestDuration] += P_df[ttau];
    
    fProdEnd = EIPDuration-mosqRestDuration;
    for (size_t n = mosqRestDuration+1; n <= fProdEnd; ++n) {
      size_t tn = (day-n)%N_v_length;
      fArray[n] =
          P_df[tn] * fArray[n - mosqRestDuration]
          + P_A[tn] * fArray[n-1];
    }
    
    
    ts = ts % N_v_length;	// index day - theta_s
    S_v[t] = P_dif[ts] * fArray[EIPDuration-mosqRestDuration] * (N_v[ts] - O_v[ts])
        + sum
        + P_A[t1]*S_v[t1]
        + P_df[ttau]*S_v[ttau];
    //END S_v
    
    
    partialEIR += S_v[t] * P_Ai_base;
  }
}


/*****************************************************************************
 *  The following code all concerns calculating the mosquito emergence rate  *
 *****************************************************************************/


/* 2008.10.06 (Nakul Chitnis)
 * We are now dealing with code that has been fully converted to C. 
 * We probably need to redo some of this work - and make sure all
 * arrays are in C format. I think we should also avoid 2-dimensional
 * arrays in C and simply use gsl_matrices so we avoid the stupid
 * C convention of the first index refering to the column.
 *
 * The code probably needs some cleaning up to remove the remains of the
 * Fortran-C interations.
 */

/* We are currently trying to create arrays and two dimensional arrays.
 * It may make more sense to just create gsl_vectors and gsl_matrices.
 * This may be a problem when we define Upsilon as a three dimensional
 * array. Maybe there is some way of dealing with that in gsl - but
 * perhaps not. Let's deal with that later.
 */ 

/* We use the naming convention that all arrays and matrices that come from 
 * Fortran and will be sent back to Fortran begin with 'F'. All vectors
 * and matrices that are created and used by C begin with 'C'.
 * Hopefully this will help to keep things less confusing
 * - although certainly not eliminate the confusion...
 */

/* In C, the first index refers to the column and the second index to the 
 * row. This is stupid, but we have no choice. However, in our attempt to
 * be as consistent as possible, we always refer to the row by 'i' and to
 * the column by 'j'. We refer to the element in the i^{th} row and j^{th} 
 * column of matrix, A, by A(j,i). This is certainly not perfect, but 
 * perhaps the best that we can do.
 */ 


void VectorTransmissionSpecies::calMosqEmergeRate (int populationSize, vector<double>& kappa) {
  /* Number of type of malaria-susceptible hosts. 
  Dimensionless.
  $m$ in model. Scalar.
  Is equal to 1 in initalization. */
  int nMalHostTypesInit = 1;
  
  /* Number of types of hosts. 
  Dimensionless.
  $n$ in model. Scalar.
  Is equal to 1 in initalization. */
  int nHostTypesInit = nMalHostTypesInit;	// TODO: plus num non-human host types

  /* Infectivity of hosts to mosquitoes.
  $K_vi$ in model. Matrix of size $n \times \theta_p$.
  In initialization, there is only one time of host.
  This is taken directly from kappa(i) from
  entomology. */
  double humanInfectivityInit[daysInYear];
  
  
  
  /**************************************************************!
  ** We enter parameters here that will later be moved to .xml **!
  ***************************************************************!
  Note that these parameters are only for one group of humans.
  We need to get parameters for multiple groups later.
  Some parameters - that should be vectors of length daysInYear
  we will enter as scalars and assume they are fixed over the 
  year. Although, in the theoretical model, we assume that they
  may vary periodically over one year, we code them as though 
  they are constant. We can change that later if they need to 
  vary. */
  

  /**************************************************************! 
  ***************************************************************!
  *********** Other main input and output parameters ************!
  ***************************************************************!*/
  
  /* The entomological inoculation rate.
  Units: 1/Time
  $\Xi_i$ in model. Matrix of size $n \times \theta_p$.
  We use this to calculate the mosquito emergence rate.
  During the initialization, it a vector with the length of the 
  annual period. */
  vector<double> EIRInit(daysInYear, 0.0);

    /*
  ***************************************************************!
  ***************************************************************!
  *************** Begin code of subroutine here *****************!
  ***************************************************************!
    */
  
  /* Now we have to deal with 
  - humanInfectivityInit - we take from kappa
  - EIRinit - we take from EIR.
  ************ Make sure that EIR is the correct one. ***********!
  I'm not sure how we should check this - but we should look into
  it more carefully at some point in time. For now, we can 
  continue with this.
  We first create arrays of length intervalsPerYear for all three.
  We then convert them to length daysInYear.
  Save the human infectivity to mosquitoes from simulation of one
  lifetime to kappa. We will then use this to create 
  humanInfectivityInit (of size daysInYear). */
  
  /* NOTE: old method (save speciesEIR from init and expand)
  // EIR should have been calculated from fourier coefic. and used for the
  // updateOneLifespan period as a forced EIR. We should now copy EIR to an
  // array of length daysInYear, without smoothing (since humans only
  // experience 1 EIR value per timestep).
  convertLengthToFullYear(EIRInit, speciesEIR);
  // The other option would be to smooth (or recalculate from FCEIR).
  //logDFTThreeModeSmooth(EIRInit, speciesEIR, daysInYear, Global::intervalsPerYear);
  
  speciesEIR.clear();	// this is finished with; frees memory?
  */
  calcInverseDFTExp(EIRInit, FCEIR);
  if(EIRRotateAngle != 0.0)
    rotateArray(EIRInit, EIRRotateAngle);
  
  
  convertLengthToFullYear(humanInfectivityInit, kappa);
  PrintArray ("kappa", humanInfectivityInit, daysInYear);
  
  /* Find an initial estimate of the mosquito emergence rate, in
  mosqEmergeRate.
  Units: Mosquitoes/Time
  
  Not explicitly included in model. Vector of length $\theta_p$.
  This initial estimate is used by a root finding algorithm
  to calculate the mosquito emergence rate. */

  /* Set values of the initial estimate for the mosquito emergence rate.
  If we have already calcuated the mosquito emergence rate for the
  given parameters separately, we can simply use that (and later
  test the resulting EIR if it matches.
  The file is assumed to contain a moquito emergence rate of length
  daysInYear.
  Otherwise, we just use a multiple of the EIR. The value of this 
  vector may not be very important, but it may speed up the root
  finding algorithm (but probably not).
  (2008.10.20: It appears to make no difference to the speed.)
    */
  ifstream file (emergenceRateFilename.c_str());
  if(file.good()) {	// file exists since it opened succesfully
    for (int i = 0; i < daysInYear; i++){
	  file >> mosqEmergeRate[i];
    }
    cout << "Read emergence rates from file: " << mosqEmergeRate[0] << ", " << mosqEmergeRate[1] << "..." << endl;
  }else{
    double temp = populationSize*populationSize*entoAvailability;
    for (int i = 0; i < daysInYear; i++) {
      mosqEmergeRate[i] = EIRInit[i]*temp;
    }
    cout << "Guessed emergence rates: " << mosqEmergeRate[0] << ", " << mosqEmergeRate[1] << "..." << endl;
  }
  file.close();
  
  // Now calculate the emergence rate.
  // The routine should finish quickly if the emergence rate is already accurate.
  if(true) {
    CalcInitMosqEmergeRate(populationSize,
                           nHostTypesInit,
                           nMalHostTypesInit,
                           humanInfectivityInit, EIRInit);
    
    // Now we've calculated the emergence rate, save it:
    ofstream file (emergenceRateFilename.c_str());
    for (int i = 0; i < daysInYear; ++i)
      file << mosqEmergeRate[i] << endl;
    file.close();
  }
}

void VectorTransmissionSpecies::convertLengthToFullYear (double FullArray[daysInYear], vector<double>& ShortArray) {
  for (size_t i=0; i < Global::intervalsPerYear; i++) {
    for (int j=0; j < Global::interval; j++) {
      FullArray[i*Global::interval + j] = ShortArray[i];
    }
  }
}



/***************************************************************************
 ************************ START SUBROUTINES HERE ***************************
 ***************************************************************************/
double VectorTransmissionSpecies::CalcInitMosqEmergeRate(int populationSize,
                                             int nHostTypesInit,
                                             int nMalHostTypesInit,
                                             double* FHumanInfectivityInitVector,
                                             vector<double>& FEIRInitVector)
{
/* Note that from here on we use the notation from "A Mathematical Model for the
  * Dynamics of Malaria in Mosquitoes Feeding on a Heterogeneous Host Population",
  * and (the publication with the periodic model - yet to be written).
  *
  * While, this may not be the easiest notation to read for someone not familiar
  * with the model, it will be easier to go directly from the equations in the paper
  * to the equations, as they will be written in the code. Since the equations are not
  * obvious in any case, anyone who wants to go through this code, will need to go 
  * through the paper, so I think that will be ok.
  *
  * There are also a number of variables defined that are difficult to describe
  * physically which we use in intermediate equations. We try to give names that
  * we use in the papers referenced above. 
  *
  *  - Any complaints about this notation (or anything else in general) can be directed 
  * to itsupport-sti@stimail.ch
  *
  * Once the paper on the periodic model is written/published - we should also include
  * the equation numbers as that may help.l
  *
  * We may append a 'CV' or 'CM' to gsl_vectors and gsl_matrices to distinguish them
  * if we feel it is necessary.
  *
  * As far as possible, we try to use gsl_vectors instead of arrays to allow more
  * flexibility. */

  // Alternative notation used:
  size_t theta_p = daysInYear;
  size_t tau = mosqRestDuration;
  size_t theta_s = EIPDuration;
  // n,m are not aliased but are: nHostTypesInit,nMalHostTypesInit

# define N_i		populationSize
# define alpha_i	entoAvailability
# define mu_vA		mosqSeekingDeathRate
# define theta_d	mosqSeekingDuration
# define P_B_i		probMosqBiting
# define P_C_i		probMosqFindRestSite
# define P_D_i		probMosqSurvivalResting
# define P_E_i		probMosqSurvivalOvipositing


  // Parameters that help to describe the order of the system.
  // Ask not why we call mt, mt. We use mt to index the system.
  // It is the maximum number of time steps we go back for $N_v$ and $O_v$.
  size_t mt = theta_s + tau -1;
  size_t eta = 2*mt + tau;	// $\eta$: The order of the system.

  
  // The set of theta_p matrices that determine the dynamics of the system
  // from one step to the next, that is, the system is described by,
  // $x(t) = \Upsilon(t) x(t-1) = \Lambda(t)$.
  // $\Upsilon(t)$ is defined over time, $1 \leq t \leq \theta_p$, 
  // where $t \in \mathbb{N}$.
  gsl_matrix** Upsilon = (gsl_matrix**) malloc(theta_p*sizeof(gsl_matrix*));

  // The set of theta_p vectors that determine the forcing of the system
  // at every time step.
  // $\Lambda(t)$ is defined over time, $1 \leq t \leq \theta_p$, 
  // where $t \in \mathbb{N}$.
  gsl_vector** Lambda = (gsl_vector**) malloc(theta_p*sizeof(gsl_vector*));

  // The periodic orbit of all eta state variables.
  gsl_vector** x_p = (gsl_vector**) malloc(theta_p*sizeof(gsl_vector*));

  // This initially contains the initial estimate of the mosquito emergence
  // rate. This is used by the root finding algorithm to calculate N_v0.
  gsl_vector* N_v0 = gsl_vector_calloc(theta_p);	// mosqEmergeRate
  memcpy (N_v0->data, mosqEmergeRate, theta_p * sizeof (*mosqEmergeRate));
  
  PrintVector("Nv0", N_v0, theta_p);

  gsl_vector* K_vi = gsl_vector_calloc(theta_p);	// humanInfectivity
  memcpy (K_vi->data, FHumanInfectivityInitVector, theta_p * sizeof (*FHumanInfectivityInitVector));
  
  // Output Parameters (for the model):
  gsl_vector* Xi_i = gsl_vector_calloc(theta_p);	// EIR
  memcpy (Xi_i->data, &FEIRInitVector[0], theta_p * sizeof (FEIRInitVector[0]));
  
  // The number of infectious mosquitoes over every day of the cycle.
  // calculated from the EIR data.
  // $S_v$ (from EIR).
  gsl_vector* S_vFromEIR = gsl_vector_calloc(theta_p);
  
  // The difference between S_vFromEIR and SvfromNv0.
  gsl_vector* S_vDiff = gsl_vector_calloc(theta_p);
  
  // State variables.
  
  // The periodic values of the total number of host-seeking mosquitoes.
  gsl_vector* Nvp = gsl_vector_calloc(theta_p);	// $N_v^{(p)}(t)$.
  // The periodic values of the number of infected host-seeking mosquitoes.
  gsl_vector* Ovp = gsl_vector_calloc(theta_p);	// $O_v^{(p)}(t)$.
  // The periodic values of the number of infectious host-seeking mosquitoes.
  gsl_vector* Svp = gsl_vector_calloc(theta_p);	// $S_v^{(p)}(t)$.
  
  // Allocate memory for gsl_matrices and initialize to 0.
  
  // $X_{\theta_p}$.
  // The product of all the evolution matrices.
  // $X_{\theta_p} = X(\theta_p+1,1)$. 
  // Refer to Cushing (1995) and the paper for the periodic entomological model
  // for more information.
  gsl_matrix* X_t_p = gsl_matrix_calloc(eta, eta);
  
  // $(\mathbb{I}-X_{\theta_p})^{-1}$.
  // The inverse of the identity matrix minus X_t_p.
  gsl_matrix* inv1Xtp = gsl_matrix_calloc(eta, eta);

# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  // We now try to print these parameters to file to make sure that 
  // they show what we want them to show.
  PrintParameters(theta_p, tau, theta_s, nHostTypesInit, nMalHostTypesInit, N_i, alpha_i,
                  mu_vA, theta_d, P_B_i, P_C_i, P_D_i, P_E_i, K_vi, Xi_i);
  // The parameter values look correct.
# endif
  
  // Derived Parameters
  
  /// Probability that a mosquito survives one day of 
  /// host-seeking but does not find a host. 
  double P_A = 0;

  /// Probability that on a given day, a mosquito finds a host of type $i$.
  /// For now, we assume that this  is  a double: 
  /// - no dependence on the phase of the period - or the type of host.
  double P_Ai = 0;	// $P_{A^i}$

  
  // Create matrices in Upsilon.
  // We also define P_A and P_Ai in the same routine. 
  // For now, we treat P_A and P_Ai as scalars since we are 
  // defining most parameters as scalars. If we do change things later, which we
  // may, then we will change the code accordingly. We will need to go through
  // a lot of changes anyway. 
  CalcUpsilonOneHost(Upsilon, &P_A, &P_Ai, theta_p, eta, mt, tau, theta_s, 
                     nHostTypesInit, nMalHostTypesInit, N_i, alpha_i, mu_vA,
                     theta_d, P_B_i, P_C_i, P_D_i, P_E_i, K_vi);


  // Calculate $X_{\theta_p}$.
  // Refer to Cushing (1995) and the paper for the periodic entomological model
  // for more information.
  FuncX(X_t_p, Upsilon, theta_p, 0, eta);

# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  char xtpname[15] = "X_t_p";
  PrintMatrix(xtpname, X_t_p, eta, eta);
# endif

  // We should now find the spectral radius of X_t_p and show that it's less than 1.
  double srXtp = CalcSpectralRadius(X_t_p, eta);

  // printf("The spectral radius of X_t_p = %e\n", srXtp);
  // getchar();

  // If the spectral radius of X_t_p is greater than or equal to 1, then
  // we are not guaranteed the existence of a unique globally asymptotically
  // stable periodic orbit; thus it does not make sense to try to match the EIR
  // for this periodic orbit.
  // 
  // For this model, all the eigenvalues should be in the unit circle. However,
  // as we cannot show that analytically, we need to check it numerically.
  if(srXtp >= 1.0){
    ostringstream msg;
    msg << "The spectral radius of X_t_p = "<< srXtp << "; expected to be less than 1.\n";
    msg << "Warning: No globally asymptotically stable periodic orbit. \n";
    msg << "Warning: All results from the entomologoical model may be meaningless. \n";
    throw xml_scenario_error (msg.str());
  }

  // Calculate the inverse of (I-X_t_p).
  CalcInv1minusA(inv1Xtp, X_t_p, eta);

# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  char inv1Xtpname[15] = "inv1minusXtp";
  PrintMatrix(inv1Xtpname, inv1Xtp, eta, eta);
# endif

  // Calculate the number of infectious host-seeking mosquitoes for the given EIR.
  CalSvfromEIRdata(S_vFromEIR, P_Ai, P_B_i, N_i, Xi_i);

# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  char SvfromEIRname[15] = "S_vFromEIR";
  PrintVector(SvfromEIRname, S_vFromEIR, theta_p);
# endif

  /**** We check the initial value of the mosquito emergence rate *****
  * If the resulting proportion of infectious host-seeking mosquitoes
  * matches that calculated from the EIR, we do not need to do any 
  * root-finding. 
  * 
  * There should probably be a clean way of doing this through the XML
  * file but for now this is probably ok.
  */
  // N_v0 already contains our estimate, but we may need a copy for root finding.
  // Input vector for root-finding:
  gsl_vector* xrootfind = gsl_vector_calloc(theta_p);
  gsl_vector_memcpy(xrootfind, N_v0);
  

  CalcSvDiff(S_vDiff, S_vFromEIR, Upsilon, N_v0, inv1Xtp, 
             eta, mt, theta_p);
# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  char InitSvDiffname[20] = "InitSvDifference";
  PrintVector(InitSvDiffname, S_vDiff, theta_p);
# endif
  double SvDiff1norm = gsl_blas_dasum(S_vDiff);	//The $l^1$ norm of S_vDiff.
  printf("The $l^1$ norm of S_vDiff is %e \n", SvDiff1norm);

  // Maximum $l^1$ distance of error of root-finding algorithm
  const double EpsAbsRF = 1.0;
  
  if(SvDiff1norm>EpsAbsRF){
    printf("The difference in Sv is greater than the tolerance. \n");
    
    if (!(Global::clOptions & CLO::ENABLE_ERC))
      throw runtime_error ("Cannot recalculate: emergence rate calculations are not enabled.");
    
    /************* We initialize variables for root-finding. **************/
    printf("Starting root-finding \n");

    // Parameters for root-finding function.
    // pararootfind = {S_vFromEIR, Upsilon, inv1Xtp, eta, mt, theta_p};
    struct SvDiffParams pararootfind;
    pararootfind.S_vFromEIR = S_vFromEIR;
    pararootfind.Upsilon = Upsilon;
    pararootfind.inv1Xtp = inv1Xtp;
    pararootfind.eta = eta;
    pararootfind.mt = mt;
    pararootfind.thetap = theta_p;

    // Set root-finding function.
    // frootfind = {&CalcSvDiff_rf, theta_p, &pararootfind};
    gsl_multiroot_function frootfind;
    frootfind.f = &CalcSvDiff_rf;
    frootfind.n = theta_p;
    frootfind.params = &pararootfind;

    // Set type of root-finding algorithm.
    const gsl_multiroot_fsolver_type* Trootfind =
        gsl_multiroot_fsolver_hybrids;
    // Allocate memory for root-finding workspace.
    gsl_multiroot_fsolver* srootfind =
        gsl_multiroot_fsolver_alloc(Trootfind, theta_p);

    printf("About to set root-finding solver \n");
    // Initialize root-finding.
    gsl_multiroot_fsolver_set(srootfind, &frootfind, xrootfind);
    printf("Set root-finding \n");

    // Print initial state (to screen):
    char fnamerootfindoutput[30] = "output_rootfinding.txt";
    PrintRootFindingStateTS(0, srootfind, theta_p, fnamerootfindoutput);
    
    // Maximum number of iterations of root-finding algorithm.
    const size_t maxiterRF = 1000;
    int status = GSL_CONTINUE;
    for (size_t iter = 1; status == GSL_CONTINUE && iter < maxiterRF; ++iter) {
      status = gsl_multiroot_fsolver_iterate(srootfind);
      PrintRootFindingStateTS(iter, srootfind, theta_p, fnamerootfindoutput);

      // Check to see if solver is stuck
      if (status){
        break;
      }

      status = gsl_multiroot_test_residual(srootfind->f, EpsAbsRF);
    }

    // Print status
    printf("status = %s \n", gsl_strerror(status)); 

    // Copy solution for N_v0 into N_v0.
    gsl_vector_memcpy(N_v0, srootfind->x);

#   ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
    char finalNv0name[15] = "FinalNv0";
    char finalSvDiffname[15] = "FinalSvDiff";
    PrintVector(finalNv0name, N_v0, theta_p);
    PrintVector(finalSvDiffname, srootfind->f, theta_p);
#   endif

		// Free memory.
    gsl_vector_free(xrootfind);
    gsl_multiroot_fsolver_free(srootfind);
  }


  // Calculate final periodic orbit.
  CalcLambda(Lambda, N_v0, eta, theta_p);
  CalcXP(x_p, Upsilon, Lambda, inv1Xtp, eta, theta_p);

  // Retrieve the periodic orbits for Nv, Ov, and Sv.
  size_t indexSv = 2*mt;
  for (size_t i=0; i<theta_p; i++){
    double temp = gsl_vector_get(x_p[i], 0);
    gsl_vector_set(Nvp, i, temp);

    temp = gsl_vector_get(x_p[i], mt);
    gsl_vector_set(Ovp, i, temp);

    temp = gsl_vector_get(x_p[i], indexSv);
    gsl_vector_set(Svp, i, temp);
  }
  
  if (Simulation::simulationTime*Global::interval < N_v_length || N_v_length > (int)theta_p)
    throw xml_scenario_error ("Initialization phase or theta_p too short");
  // For day over N_v_length days prior to the next timestep's day.
  int endDay = (Simulation::simulationTime+1) * Global::interval;
  for (int day = endDay - N_v_length; day < endDay; ++day) {
    size_t t = day % N_v_length;
    size_t i = (theta_p + day - endDay) % theta_p;
    
    N_v[t] = gsl_vector_get (Nvp, i);
    O_v[t] = gsl_vector_get (Ovp, i);
    S_v[t] = gsl_vector_get (Svp, i);
  }
  
  
# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  char Nvpname[15] = "NvPO";
  char Ovpname[15] = "OvPO";
  char Svpname[15] = "SvPO";
  PrintVector(Nvpname, Nvp, theta_p);
  PrintVector(Ovpname, Ovp, theta_p);
  PrintVector(Svpname, Svp, theta_p);
# endif

  for (size_t i=0; i<theta_p; i++){
    gsl_vector_free(Lambda[i]);
    gsl_vector_free(x_p[i]);
  }


  // Copy the mosquito emergence rate to the C array.
  memcpy (mosqEmergeRate, N_v0->data, theta_p * sizeof (*mosqEmergeRate));

  // Deallocate memory for vectors and matrices.
  for (size_t i=0; i<theta_p; i++){
    gsl_matrix_free(Upsilon[i]);
  }

  gsl_vector_free(N_v0);
  gsl_vector_free(K_vi);
  gsl_vector_free(Xi_i);
  gsl_vector_free(S_vFromEIR);
  gsl_vector_free(S_vDiff);
  gsl_vector_free(Nvp);
  gsl_vector_free(Ovp);
  gsl_vector_free(Svp);
	

  gsl_matrix_free(X_t_p);
  gsl_matrix_free(inv1Xtp);

  free(Upsilon);
  free(Lambda);
  free(x_p);

  return 0.0;
}

void VectorTransmissionSpecies::logDFTThreeModeSmooth (double* smoothArray,
                                               double* originalArray,
                                               int SALength,
                                               int OALength)
{
  // Frequency
  double foa = 1.0/OALength;
  double woa = 2*M_PI * foa;
  double wsa = 2*M_PI/SALength;
  
  double tempsuma0 = 0.0;
  double tempsuma1 = 0.0;
  double tempsumb1 = 0.0;
  double tempsuma2 = 0.0;
  double tempsumb2 = 0.0;
  
  // Calculate first three Fourier modes
  for (int t=0; t < OALength; t++) {
    double yt = log(originalArray[t]);
    double woa_t = woa*t;
    tempsuma0 += yt;
    tempsuma1 += (yt*cos(woa_t));
    tempsumb1 += (yt*sin(woa_t));
    tempsuma2 += (yt*cos(2*woa_t));
    tempsumb2 += (yt*sin(2*woa_t));       
  }
  
  // Fourier Coefficients
  double a0 = (  foa)*tempsuma0;
  double a1 = (2*foa)*tempsuma1;
  double b1 = (2*foa)*tempsumb1;
  double a2 = (2*foa)*tempsuma2;
  double b2 = (2*foa)*tempsumb2;
  
  // Calculate inverse discrete Fourier transform
  for (int t=0; t < SALength; t++){
    double wsa_t = wsa*(t+1);
    smoothArray[t] = 
        exp(a0 + a1*cos(wsa_t) + b1*sin(wsa_t) + 
        a2*cos(2*wsa_t) + b2*sin(2*wsa_t));
  }

#   ifdef TransmissionModel_PrintSmoothArray
  PrintArray("SmoothArray", smoothArray, SALength);
#   endif
}



void VectorTransmissionSpecies::calcInverseDFTExp(vector<double>& tArray, vector<double>& FC) {
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

void VectorTransmissionSpecies::rotateArray(vector<double>& rArray, double rAngle) {
  vector<double> tempArray (rArray.size());
  size_t rotIndex = size_t ((rAngle*rArray.size())/(2.0*M_PI));

  for (size_t i=0; i < rArray.size(); i++) {
    tempArray[(i+rotIndex) % rArray.size()] = rArray[i];
  }
  
# ifdef TransmissionModel_PrintRotateArray
  PrintArray((char*)"PrerotationArray", rArray);
  PrintArray((char*)"PostrotationArray", tempArray);
# endif
  rArray = tempArray;
}
