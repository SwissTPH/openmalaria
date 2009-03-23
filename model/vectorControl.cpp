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

// Variable names largely come from Nakul Chitnis's paper:
// "A mathematical model for the dynamics of malaria in mosquitoes feeding on
// a heterogeneous host population" (3rd Oct. 2007).

/* Entomology model coordinator: Nakul Chitnis. */ 

#include "vectorControl.h"
#include "VectorControlInternal.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>


/*****************************************************************************
 *                       New code, written in C++                            *
 *****************************************************************************/


VectorControl::VectorControl () {
  //TODO: initialize these and perhaps other params properly:
  // Get params from xml:
  //alpha = ;
  /* This will be defined in the xml file as a real.
  We will define it as a integer - and convert after multiplying
  by the length of the Global::interval. Scalar. */
  mosqRestDuration=3;		// TODO: Move to XML

  mosqSeekingDeathRate=1.6;	// TODO: Move to XML
  mosqSeekingDuration=0.33;	// TODO: Move to XML
  
  probMosqEggLaying=0;
  
  if (EIPDuration < 1 || mosqRestDuration < 1) {
    cerr << "Code expects EIPDuration and mosqRestDuration to be a positive number of days" << endl;
    throw 0;
  }
  
  N_v_length = EIPDuration + mosqRestDuration;
  
  P_A = new double[N_v_length * 6];	// allocate memory for all arrays
  size_t l = N_v_length*sizeof(double);
  P_df	= P_A + l;
  P_dif	= P_df + l;
  N_v	= P_dif + l;
  O_v	= N_v + l;
  S_v	= O_v + l;
  
  //TODO: initialise arrays (run initial simulation for N_v_length-1 days?)
}
VectorControl::~VectorControl () {
  delete[] P_A;
}

void VectorControl::initMainSimulation(int populationSize) {
  cerr << "Warning: using incomplete VectorControl transmission model!" << endl;
  calMosqEmergeRate (populationSize);
}

// dummy functions until they're implemented:
double VectorControl::getExpectedNumberOfInfections (Human& human, double age_adj_EIR) {
  // I'm not really sure what this should do (hardy).
  cerr << "dummy function getExpectedNumberOfInfections called" << endl;
  return 0.0;
}

/** Calculate EIR for host, using the fixed point of difference eqns. */
double VectorControl::calculateEIR(int simulationTime, Human& host) {
  /* Calculates EIR per individual (hence N_i == 1).
   *
   * See comment in advancePeriod for method. */
  
  return partialEIR
      * host.entoAvailability()
      * host.probMosqSurvivalBiting();	// probability of biting, once commited
}


// Per 5-day step:
void VectorControl::advancePeriod (const std::list<Human>& population, int simulationTime) {
  /* Largely equations correspond to Nakul Chitnis's model in
    "A mathematic model for the dynamics of malaria in
    mosquitoes feeding on a heterogeneous host population"
  section 2, 3.5-3.6, plus extensions to a non-autonomous case (where
  emergence rate varies over the year).
  
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
  
  // Per Global::interval (hosts don't update per day):
  double totalAvailability = 0.0;
  for (std::list<Human>::const_iterator h = population.begin(); h != population.end(); ++h)
    totalAvailability += h->entoAvailability();
  
  
  // Summed per day:
  partialEIR = 0.0;
  
  
  // The code within the for loop needs to run per-day, wheras the main
  // simulation uses Global::interval day (currently 5 day) time steps.
  int endDay = (simulationTime+1) * Global::interval;
  for (int day = simulationTime * Global::interval; day < endDay; ++day) {
    // Indecies for today, yesterday and mosqRestDuration days back:
    size_t t    = day % N_v_length;
    size_t t1   = (t - 1) % N_v_length;
    size_t ttau = (t - mosqRestDuration) % N_v_length;
    
    double leaveHostRate = totalAvailability + mosqSeekingDeathRate;
  
    // Probability of a mosquito not finding a host this day:
    P_A[t] = exp(-leaveHostRate * mosqSeekingDuration);
    
    double P_Ai_base = (1.0 - P_A[t]) / leaveHostRate;
      // NOTE: already have kappa array - is same?
    
    // NC's non-autonomous model provides two methods for calculating P_df and
    // P_dif; here we assume that P_E is constant.
    double sum = 0.0;
    double sum_dif = 0.0;
    for (std::list<Human>::const_iterator h = population.begin(); h != population.end(); ++h) {
      double prod = h->entoAvailability() * h->probMosqSurvivalBiting() * h->probMosqSurvivalResting();
      sum += prod;
      sum_dif += prod;//FIXME: * h->K_vi();
    }
    P_df[t] = sum * P_Ai_base * probMosqEggLaying;
    P_dif[t] = sum_dif * P_Ai_base * probMosqEggLaying;
    
    
    //FIXME: formulas need adjusting for NC's non-autonomous model
    N_v[t] = mosqEmergeRate[day%daysInYear]
        + P_A[t]  * N_v[t1]
        + P_df[t] * N_v[ttau];
    O_v[t] = P_dif[t] * (N_v[ttau] - O_v[ttau])
        + P_A[t]  * O_v[t1]
        + P_df[t] * O_v[ttau];
    
    sum = 0.0;	// Sums making S_v[t] in eqn. (3c)
    int k_p = EIPDuration/mosqRestDuration - 1;		// k_+
    for (int j = 0; j <= k_p; ++j) {
      int temp = EIPDuration-(j+1)*mosqRestDuration;
      sum += gsl_sf_choose(temp+j, j)
          * gsl_pow_int(P_A[t] , temp)
          * gsl_pow_int(P_df[t], j);
    }
    size_t ts = (t - EIPDuration) % N_v_length;
    S_v[t] = P_dif[t] * sum * N_v[ts] - O_v[ts]
        + P_A[t]  * S_v[t1]
        + P_df[t] * S_v[ttau];	// + second sum:
    
    sum = 0.0;
    for (int l = 1; l < mosqRestDuration; ++l) {
      double S_subsum = 0.0;
      k_p = (EIPDuration+l)/mosqRestDuration - 2;	// k_l+
      for (int j = 0; j <= k_p; ++j) {
        int temp = EIPDuration+l-(j+2)*mosqRestDuration;
        S_subsum += gsl_sf_choose(temp+j, j)
            * gsl_pow_int(P_A[t] , temp)
            * gsl_pow_int(P_df[t], j);
      }
      ts = (t - EIPDuration - l) % N_v_length;
      sum += S_subsum * (N_v[ts] - O_v[ts]);
    }
    S_v[t] += sum * P_df[t];
    
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


void VectorControl::calMosqEmergeRate (int populationSize) {
  /* Number of types of hosts. 
  Dimensionless.
  $n$ in model. Scalar.
  Is equal to 1 in initalization. */
  int nHostTypesInit = 1;

  /* Number of type of malaria-susceptible hosts. 
  Dimensionless.
  $m$ in model. Scalar.
  Is equal to 1 in initalization. */
  int nMalHostTypesInit = 1;
  
  /* Infectivity of hosts to mosquitoes.
  $K_vi$ in model. Matrix of size $n \times \theta_p$.
  In initialization, there is only one time of host.
  This is taken directly from initialKappa(i) from
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
  vary.
  
  We assume that we will read them from the xml file as doubles
  so that they will be passed into Fortran as real*8. We 
  therefore define them here as real*8. The results, then, when
  used in the main simulation can be converted to reals.
  When they are defined in the xml file, we can rearrange 
  parameters here so that one group is taken directly from the
  xml file and the other group from other Fortran routines.
  And we still initialize the parameters here, but do not
  define them. */
  
  /* Availability rate of hosts to mosquitoes.
  Units: 1/(Animals * Time)
  $\alpha_i$ in model. Matrix of size $n \times \theta_p$.
  In initialization, there is only one time of host.
  We also assume that this does not change over the cycle.
  This adds a further degree of freedom - the value should be 
  chosen carefully as there will be an inverse relationhip
  between $\sum_{i=1}^n \alpha_i N_i$ and $N_{v0}$. */
  // We set the host availability relative to the population size.
  double hostAvailabilityRateInit=7.0/populationSize;	// TODO: Move absolute availability to XML
    
  /*! Probability of a mosquito biting a host given that it has encountered the host.
  
  Dimensionless.
  $P_{B_i}$ in model. Matrix of size $n \times \theta_p$.
  For now, we assume this does not change over the cycle. */
  double mosqProbBiting=0.95;		// TODO: Move to XML

  /*! Probability of a mosquito finding a resting side given that it has bitten a host.
  
  Dimensionless.
  $P_{B_i}$ in model. Matrix of size $n \times \theta_p$.
  For now, we assume this does not change over the cycle. */
    
  double mosqProbFindRestSite=0.95;	// TODO: Move to XML
    
  /*! Probability of a mosquito surviving the resting period given that it has found a resting site.
   * 
  Dimensionless.
  $P_{D_i}$ in model. Matrix of size $n \times \theta_p$.
  For now, we assume this does not change over the cycle. */
  double mosqProbResting=0.94;		// TODO: Move to XML
  
  /*! Probability of a mosquito ovipositing and returning to host-seeking given that it has survived the resting period.
   * 
  Dimensionless.
  $P_{E_i}$ in model. Matrix of size $n \times \theta_p$.
  For now, we assume this does not change over the cycle. */
  double mosqProbOvipositing=0.93;	// TODO: Move to XML

  int ifUseNv0Guess = 0;	// TODO: Move to XML.
				// Use a predefined array for the initial
				// mosquito emergence rate.
				// Perhaps calculated in a different
				// iteration.
  char Nv0guessfilename[20] = "N_v0-Initial.txt";
				// TODO: Move to XML.
				// File that contains our initial guess for the 
				// mosquito emergence rate.
				// To be used only if ifUseNv0Guess=1.

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
  double EIRInit[daysInYear];

    /*
  ***************************************************************!
  ***************************************************************!
  *************** Begin code of subroutine here *****************!
  ***************************************************************!
    */
  
  /* Now we have to deal with 
  - humanInfectivityInit - we take from initialKappa
  - EIRinit - we take from EIR.
  ************ Make sure that EIR is the correct one. ***********!
  I'm not sure how we should check this - but we should look into
  it more carefully at some point in time. For now, we can 
  continue with this.
  We first create arrays of length intervalsperYear for all three.
  We then convert them to length daysInYear.
  Save the human infectivity to mosquitoes from simulation of one
  lifetime to initialKappa. We will then use this to create 
  humanInfectivityInit (of size daysInYear). */

	// We need to decide how we deal with the EIR - if we smooth it out
	// over the entire year, or leave it constant over the Global::interval length.
	// I think it would be better ot smooth it over the full year. 
	// It may not be fully accurate - but we are in any case going to lose
	// some accuracy over the difference between the time step of the human
	// simulation model and the mosquito transmission model.
	// Note that we smooth over the entire year, we are slightly shifting
	// the EIR a little bit to the right.
	// We have a lot of if statements and flags here. We need to clean this 
	// up eventually when the flags are moved to XML.
  if( true ){
    
    if(ifUseFC){
      printf("Calculating inverse discrete Fourier transform over full year \n");
      calcInverseDFTExp(EIRInit, daysInYear, FCEIR, FCEIRX);
    } else if(FTSmoothEIR==1){
      printf("Smoothing and expanding EIR \n");
      logDFTThreeModeSmooth(EIRInit, origEIR, daysInYear, Global::intervalsPerYear);
    } else
      convertLengthToFullYear(EIRInit, EIR);
      
  } else
    convertLengthToFullYear(EIRInit, EIR);
    
  if(ifrotateEIR){
    rotateArray(EIRInit, daysInYear, EIRRotateAngle);
    printf("Rotating EIR \n");
  }

# ifdef VectorControl_PRINT_calMosqEmergeRate
  char origEIRname[15] = "OrigEIR";
  char shortEIRname[15] = "ShortEIR";
  char longEIRname[15] = "LongEIR";
  PrintArray(fnametestentopar, origEIRname, origEIR, Global::intervalsPerYear);
  PrintArray(fnametestentopar, shortEIRname, EIR, Global::intervalsPerYear);
  PrintArray(fnametestentopar, longEIRname, EIRInit, daysInYear);
# endif

  convertLengthToFullYear(humanInfectivityInit, initialKappa);
# ifdef VectorControl_PRINT_calMosqEmergeRate
  char shortKviname[15] = "ShortKvi";
  char longKviname[15] = "LongKvi";
  PrintArray(fnametestentopar, shortKviname, initialKappa, Global::intervalsPerYear);
  PrintArray(fnametestentopar, longKviname, humanInfectivityInit, daysInYear);
# endif
  
  
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
  if(ifUseNv0Guess){
    // Read from file.
    FILE* fNv0 = fopen(Nv0guessfilename, "r");
    for(int i = 0; i < daysInYear; i++){
      fscanf(fNv0, "%le", &mosqEmergeRate[i]);
    }
    fclose(fNv0);
  }else{
    double temp = populationSize*populationSize*hostAvailabilityRateInit;
    for (int i = 0; i < daysInYear; i++) {
      mosqEmergeRate[i] = EIRInit[i]*temp;
    }
  }

  // Now calculate the emergence rate:
  if(ifCalcMosqEmergeRate)
    CalcInitMosqEmergeRate(populationSize, EIPDuration,
                           nHostTypesInit,
                           nMalHostTypesInit, hostAvailabilityRateInit,
                           mosqProbBiting, mosqProbFindRestSite,
                           mosqProbResting, mosqProbOvipositing,
                           humanInfectivityInit, EIRInit);
}

void VectorControl::convertLengthToFullYear (double FullArray[daysInYear], double* ShortArray) {
  if (daysInYear != Global::interval*Global::intervalsPerYear)
    throw 0;
  
  for (int i=0; i < Global::intervalsPerYear; i++) {
    for (int j=0; j < Global::interval; j++) {
      FullArray[i*Global::interval+j] = ShortArray[i];
    }
  }
}



/***************************************************************************
 ************************ START SUBROUTINES HERE ***************************
 ***************************************************************************/
double VectorControl::CalcInitMosqEmergeRate(int populationSize,
                                             int EIPDuration,
                                             int nHostTypesInit,
                                             int nMalHostTypesInit,
                                             double hostAvailabilityRateInit,
                                             double mosqProbBiting,
                                             double mosqProbFindRestSite,
                                             double mosqProbResting,
                                             double mosqProbOvipositing,
                                             double* FHumanInfectivityInitVector,
                                             double* FEIRInitVector)
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
# define alpha_i	hostAvailabilityRateInit
# define mu_vA		mosqSeekingDeathRate
# define theta_d	mosqSeekingDuration
# define P_B_i		mosqProbBiting
# define P_C_i		mosqProbFindRestSite
# define P_D_i		mosqProbResting
# define P_E_i		mosqProbOvipositing


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
  
  gsl_vector* K_vi = gsl_vector_calloc(theta_p);	// humanInfectivity
  memcpy (K_vi->data, FHumanInfectivityInitVector, theta_p * sizeof (*FHumanInfectivityInitVector));
  
  // Output Parameters (for the model):
  gsl_vector* Xi_i = gsl_vector_calloc(theta_p);	// EIR
  memcpy (Xi_i->data, FEIRInitVector, theta_p * sizeof (*FEIRInitVector));
  
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

# ifdef VectorControl_PRINT_CalcInitMosqEmergeRate
  // We now try to print these parameters to file to make sure that 
  // they show what we want them to show.
  PrintParameters(fnametestentopar, theta_p, tau, theta_s, nHostTypesInit, nMalHostTypesInit, N_i, alpha_i,
                  mu_vA, theta_d, P_B_i, P_C_i, P_D_i, P_E_i, K_vi, Xi_i, Nv0guess);
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
                     theta_d, P_B_i, P_C_i, P_D_i, P_E_i, K_vi,
                     fnametestentopar);


  // Calculate $X_{\theta_p}$.
  // Refer to Cushing (1995) and the paper for the periodic entomological model
  // for more information.
  FuncX(X_t_p, Upsilon, theta_p, 0, eta);

# ifdef VectorControl_PRINT_CalcInitMosqEmergeRate
  char xtpname[15] = "X_t_p";
  PrintMatrix(fnametestentopar, xtpname, X_t_p, eta, eta);
# endif

  // We should now find the spectral radius of X_t_p and show that it's less than 1.
  double srXtp = CalcSpectralRadius(X_t_p, eta, fnametestentopar);

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
    printf("Warning: The spectral radius of $X_{tp}$ is not less than 1.\n"); 
    printf("Warning: No globally asymptotically stable periodic orbit. \n");
    printf("Warning: All results from the entomologoical model may be meaningless. \n");
    printf("The spectral radius of X_t_p = %e\n", srXtp);
    throw 0;
  }

  // Calculate the inverse of (I-X_t_p).
  CalcInv1minusA(inv1Xtp, X_t_p, eta, fnametestentopar);

# ifdef VectorControl_PRINT_CalcInitMosqEmergeRate
  char inv1Xtpname[15] = "inv1minusXtp";
  PrintMatrix(fnametestentopar, inv1Xtpname, inv1Xtp, eta, eta);
# endif

  // Calculate the number of infectious host-seeking mosquitoes for the given EIR.
  CalSvfromEIRdata(S_vFromEIR, P_Ai, P_B_i, N_i, Xi_i);

# ifdef VectorControl_PRINT_CalcInitMosqEmergeRate
  char SvfromEIRname[15] = "S_vFromEIR";
  PrintVector(fnametestentopar, SvfromEIRname, S_vFromEIR, theta_p);
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
             eta, mt, theta_p, fnametestentopar);
# ifdef VectorControl_PRINT_CalcInitMosqEmergeRate
  char InitSvDiffname[20] = "InitSvDifference";
  PrintVector(fnametestentopar, InitSvDiffname, S_vDiff, theta_p);
# endif
  double SvDiff1norm = gsl_blas_dasum(S_vDiff);	//The $l^1$ norm of S_vDiff.
  printf("The $l^1$ norm of S_vDiff is %e \n", SvDiff1norm);

  // Maximum $l^1$ distance of error of root-finding algorithm
  const double EpsAbsRF = 1.0;
  
  if(SvDiff1norm>EpsAbsRF){
    printf("The difference in Sv is greater than the tolerance. \n");

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

#   ifdef VectorControl_PRINT_CalcInitMosqEmergeRate
    char finalNv0name[15] = "FinalNv0";
    char finalSvDiffname[15] = "FinalSvDiff";
    PrintVector(fnametestentopar, finalNv0name, N_v0, theta_p);
    PrintVector(fnametestentopar, finalSvDiffname, srootfind->f, theta_p);
#   endif

		// Free memory.
    gsl_vector_free(xrootfind);
    gsl_multiroot_fsolver_free(srootfind);
  }


# ifdef VectorControl_PRINT_CalcInitMosqEmergeRate
  // Calculate final periodic orbit and print out values.

  // Calculate final periodic orbit.
  CalcLambda(Lambda, N_v0, eta, theta_p, fnametestentopar);
  CalcXP(x_p, Upsilon, Lambda, inv1Xtp, eta, theta_p, fnametestentopar);

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
  char Nvpname[15] = "NvPO";
  char Ovpname[15] = "OvPO";
  char Svpname[15] = "SvPO";
  PrintVector(fnametestentopar, Nvpname, Nvp, theta_p);
  PrintVector(fnametestentopar, Ovpname, Ovp, theta_p);
  PrintVector(fnametestentopar, Svpname, Svp, theta_p);

  for (size_t i=0; i<theta_p; i++){
    gsl_vector_free(Lambda[i]);
    gsl_vector_free(x_p[i]);
  }
# endif


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
