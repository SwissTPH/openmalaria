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
void VectorControl::advancePeriod (int simulationTime) {
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
  double totalAvailability = 0;// FIXME: =popSize * alpha*sum_host(host.relativeAvailability);
  double leaveHostRate = totalAvailability + mosqSeekingDeathRate;
  
  // Summed per day:
  partialEIR = 0.0;
  
  
  // The code within the for loop needs to run per-day, wheras the main
  // simulation uses Global::interval day (currently 5 day) time steps.
  int endDay = (simulationTime+1) * Global::interval;
  for (int day = simulationTime * Global::interval; day < endDay; ++day) {
    // Indecies for today, yesterday and mosqRestDuration days back:
    int t    = day % N_v_length;
    int t1   = (t - 1) % N_v_length;
    int ttau = (t - mosqRestDuration) % N_v_length;
    
    
    // Probability of a mosquito not finding a host this day:
    P_A[t] = exp(-leaveHostRate * mosqSeekingDuration);
    
    double P_Ai_base = (1.0 - P_A[t]) / leaveHostRate;
      // NOTE: already have kappa array - is same?
    
    // NC's non-autonomous model provides two methods for calculating P_df;
    // here we assume P_E is constant.
    P_df[t] =0//FIXME: =sum_host (host.entoAvailability * host.probMosqSurvivalBiting * host.probMosqSurvivalResting)
     * P_Ai_base * probMosqEggLaying;
    P_dif[t] =0//FIXME: =sum_host (host.entoAvailability * host.probMosqSurvivalBiting * host.probMosqSurvivalResting * host.K_vi)
     * P_Ai_base * probMosqEggLaying;
    
    //FIXME: formulas need adjusting for NC's non-autonomous model
    N_v[t] = mosqEmergeRate[day%daysInYear]
        + P_A[t]  * N_v[t1]
        + P_df[t] * N_v[ttau];
    O_v[t] = P_dif[t] * (N_v[ttau] - O_v[ttau])
        + P_A[t]  * O_v[t1]
        + P_df[t] * O_v[ttau];
    
    double S_sum = 0.0;	// Sums making S_v[t] in eqn. (3c)
    int k_p = EIPDuration/mosqRestDuration - 1;		// k_+
    for (int j = 0; j <= k_p; ++j) {
      int temp = EIPDuration-(j+1)*mosqRestDuration;
      S_sum += gsl_sf_choose(temp+j, j)
          * gsl_pow_int(P_A[t] , temp)
          * gsl_pow_int(P_df[t], j);
    }
    int ts = (t - EIPDuration) % N_v_length;
    S_v[t] = P_dif[t] * S_sum * N_v[ts] - O_v[ts]
        + P_A[t]  * S_v[t1]
        + P_df[t] * S_v[ttau];	// + second sum:
    
    S_sum = 0.0;
    for (size_t l = 1; l < mosqRestDuration; ++l) {
      double S_subsum = 0.0;
      k_p = (EIPDuration+l)/mosqRestDuration - 2;	// k_l+
      for (int j = 0; j <= k_p; ++j) {
        int temp = EIPDuration+l-(j+2)*mosqRestDuration;
        S_subsum += gsl_sf_choose(temp+j, j)
            * gsl_pow_int(P_A[t] , temp)
            * gsl_pow_int(P_df[t], j);
      }
      ts = (t - EIPDuration - l) % N_v_length;
      S_sum += S_subsum * (N_v[ts] - O_v[ts]);
    }
    S_v[t] += S_sum * P_df[t];
    
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
    for(size_t i=0; i<daysInYear; i++){
      fscanf(fNv0, "%le", &mosqEmergeRate[i]);
    }
    fclose(fNv0);
  }else{
    double temp = populationSize*populationSize*hostAvailabilityRateInit;
    for (size_t i=0; i<daysInYear; i++) {
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
  
  for (size_t i=0; i < Global::intervalsPerYear; i++) {
    for (size_t j=0; j < Global::interval; j++) {
      FullArray[i*Global::interval+j] = ShortArray[i];
    }
  }
}


/***************************************************************************
 *********************** STRUCTURE DEFINITIONS *****************************
 ***************************************************************************/

// Structure that contains the parameters for the function used in the 
// root-finding algorithm to find the emergence rate that matches the 
// number of infectious host-seeking mosquitoes.
struct SvDiffParams
{
  gsl_vector* S_vFromEIR;
  gsl_matrix** Upsilon;
  gsl_matrix* inv1Xtp;
  int eta;
  int mt;
  int thetap;
};



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
# define theta_p	daysInYear
# define tau		mosqRestDuration
# define theta_s	EIPDuration
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
  int mt = theta_s + tau -1;
  int eta = 2*mt + tau;	// $\eta$: The order of the system.
  int indexNv = 0;	// Index of the total number of host-seeking mosquitoes.
  int indexOv = mt;	// Index of the infected host-seeking mosquitoes.
  int indexSv = 2*mt;	// Index of the infectious host-seeking mosquitoes.

  
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
  for (size_t i=0; i<theta_p; i++){
    double temp = gsl_vector_get(x_p[i], indexNv);
    gsl_vector_set(Nvp, i, temp);

    temp = gsl_vector_get(x_p[i], indexOv);
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


void CalcUpsilonOneHost(gsl_matrix** Upsilon, double* PAPtr, 
                        double* PAiPtr, int theta_p, int eta, int mt, int tau, 
                        int theta_s, int n, int m, double N_i, double alpha_i, 
                        double mu_vA, double theta_d, double P_B_i, double P_C_i, double P_D_i, 
                        double P_E_i, gsl_vector* K_vi, char fntestentopar[])
{
	// $P_{dif}$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle and gets infected.
  gsl_vector* Pdif;
	// $P_{duf}$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle and does not get infected.
  gsl_vector* Pduf; 
	
	// The probability of a mosquito surviving the extrinsic incubation
	// period. 
	// This is the sum of j from 0 to k_+ in (2.3c).
  double sumkplus; 
  double* sumkplusPtr;

	// This is an array of the sums from 0 to k_{l+} in (2.3c).
	// Note that sumklplus here is defined as sumlv in MATLAB.
  double* sumklplus;
  sumklplus = (double *)malloc((tau-1)*sizeof(double));


	// Initialize gsl_vectors.
  Pdif = gsl_vector_calloc(theta_p);
  Pduf = gsl_vector_calloc(theta_p);


	// We note again, that this code is written assuming there is only
	// one type of hosts. 
	
	// Refer to papers noted above for equations.
        // P_A and P_Ai are described in CalcInitMosqEmergeRate.
  double P_A = exp(-(alpha_i*N_i+mu_vA)*theta_d);
  double P_Ai = (1-P_A)*(alpha_i*N_i)/(alpha_i*N_i+mu_vA);
	// $P_{df}$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle.
  double Pdf = P_Ai*P_B_i*P_C_i*P_D_i*P_E_i;

	// Evaluate Pdif and Pduf.
	// Note that these formulae are invalid for n>1.
	// We can generalize these to any n later 
	// - perhaps in a different function.
	
	// Pdif:
  gsl_vector_memcpy(Pdif, K_vi);
  gsl_vector_scale(Pdif, Pdf);

	// Pduf:
  gsl_vector_set_all(Pduf, 1.0);
  gsl_vector_sub(Pduf, K_vi);
  gsl_vector_scale(Pduf, Pdf);
	

  sumkplus = 0;
  sumkplusPtr = &sumkplus;

	// Calculate probabilities of mosquito surviving the extrinsic
	// incubation period.]
	// These currently do not depend on the phase of the period.
  CalcPSTS(sumkplusPtr, sumklplus, theta_s, tau, P_A, Pdf);
  sumkplus = *sumkplusPtr;

	/* // Print sumkplus and sumklplus
  printf("In CalcUpsilonOneHost, sumkplus = %f\n", sumkplus);
  for(l=1; l<=tau-1; l++){
  printf("sumklplus(%d) = %f \n", l, sumklplus[l-1]);
}
  getchar();
        */ 
        

	// We start creating the matrices now.
	// Refer to Section 2.1 of JBD Paper for how this matrix is created.
  for (size_t k=0; k < theta_p; k++){
    Upsilon[k] = gsl_matrix_calloc(eta, eta);

    for (size_t i=0; i<eta; i++){
			// Set 1's along the subdiagonal of all rows except the three
			// rows for the the main system variables.
      if(!((i==0) || (i==mt) || (i==(2*mt)))){
        gsl_matrix_set(Upsilon[k],i,i-1,1.0);
      }
    }

		
		// for $N_v$.
    gsl_matrix_set(Upsilon[k],0,0,P_A);
    double temp = Pdf + gsl_matrix_get(Upsilon[k], 0, tau-1);
    gsl_matrix_set(Upsilon[k],0,tau-1,temp);

		// for $O_v$.
		// We add theta_p to i, to ensure that it's positive.
		// % is the mod function.
    temp = gsl_vector_get(Pdif,(k+theta_p-tau)%theta_p);		
    gsl_matrix_set(Upsilon[k],mt,tau-1,temp);
    gsl_matrix_set(Upsilon[k],mt,mt,P_A);
    temp = gsl_vector_get(Pduf,(k+theta_p-tau)%theta_p) 
        + gsl_matrix_get(Upsilon[k], mt, mt+tau-1);	
    gsl_matrix_set(Upsilon[k],mt,mt+tau-1,temp);

		// for $S_v$.
    temp = gsl_vector_get(Pdif,(k+theta_p-theta_s)%theta_p)*sumkplus;
    gsl_matrix_set(Upsilon[k],2*mt,theta_s-1,temp);
    gsl_matrix_set(Upsilon[k],2*mt,mt+theta_s-1,-temp);
    for (size_t l=1; l <= tau-1; l++){
      temp = gsl_vector_get(Pdif,(k+theta_p-theta_s-l)%theta_p)*sumklplus[l-1];
      gsl_matrix_set(Upsilon[k],2*mt, theta_s+l-1, temp);
      gsl_matrix_set(Upsilon[k],2*mt, mt+theta_s+l-1, -temp);
    }
    gsl_matrix_set(Upsilon[k], 2*mt, 2*mt, P_A);
    temp = Pdf + gsl_matrix_get(Upsilon[k], 2*mt, 2*mt+tau-1);
    gsl_matrix_set(Upsilon[k], 2*mt, 2*mt+tau-1, temp);


  }

# ifdef VectorControl_PRINT_CalcUpsilonOneHost
  // We should try to print some of these matrices out to see what they  look like.
  PrintUpsilon(fntestentopar, Upsilon, theta_p, eta, P_A, P_Ai, Pdf, Pdif, Pduf);
# endif

	// Reference pointers.
  *PAPtr = P_A;
  *PAiPtr = P_Ai;
	
	// Deallocate memory for vectors
  gsl_vector_free(Pdif);
  gsl_vector_free(Pduf);
}


int CalcSvDiff_rf(const gsl_vector* x, void* p, gsl_vector* f){
	// Add a static variable to keep track of how often we are in this routine.
  static int counterSvDiff = 0;

	// The $l^1$ norm of Svdiff.
  double SvDiff1norm;

	// Cast the incoming pointer, p, to point at a structure of type
	// SvDiffParams.
  SvDiffParams* params = (SvDiffParams *) p;

	// Assign parameters from params to variables defined in this routine.
  gsl_vector* S_vFromEIR = (params->S_vFromEIR);
  gsl_matrix** Upsilon = (params->Upsilon);
  gsl_matrix* inv1Xtp = (params->inv1Xtp);
  int eta = (params->eta);
  int mt = (params->mt);
  int theta_p = (params->thetap);

	// It would be cleaner to read in the name of this file as an input
	// parameter but for now, we leave it out of the root-finding
	// algorithm and simply redefine it here.
  char fnametestentopar[30] = "output_ento_para.txt";	

	// Recreate a new N_v0 so that we're not restricted by const problems.
  gsl_vector* N_v0 = gsl_vector_calloc(theta_p);
  gsl_vector_memcpy(N_v0, x);


  counterSvDiff++;
  printf("In CalcSvDiff_rf for the %d th time \n", counterSvDiff);

	// To set f, we simply call CalcSvDiff. It's probably easier than rewriting
	// this code.
  CalcSvDiff(f, S_vFromEIR, Upsilon, N_v0, inv1Xtp, 
             eta, mt, theta_p, fnametestentopar);

	// Calculate the l^1 norm of f.
  SvDiff1norm = gsl_blas_dasum(f);
	
	// gsl_vector_set_all(f, 2.3);

  printf("The $l^1$ norm of S_vDiff is %e \n", SvDiff1norm);
	
  gsl_vector_free(N_v0);
  return GSL_SUCCESS;
}


void CalcSvDiff(gsl_vector* S_vDiff, gsl_vector* S_vFromEIR, 
                gsl_matrix** Upsilon, gsl_vector* N_v0, gsl_matrix* inv1Xtp, 
                int eta, int mt, int theta_p, char fntestentopar[])
{
  char SvfromNv0name[15] = "SvfromNv0";

	// The set of theta_p vectors that determine the forcing of the system
	// at every time step.
	// $\Lambda(t)$ is defined over time, $1 \leq t \leq \theta_p$, 
	// where $t \in \mathbb{N}$.
  gsl_vector** Lambda = (gsl_vector**) malloc(theta_p*sizeof(gsl_vector*));

	// The full periodic orbit.
  gsl_vector** x_p = (gsl_vector**) malloc(theta_p*sizeof(gsl_vector*));

	// Periodic orbit of the number of infectious mosquitoes calculated for
	// the given N_v0.
	// $S_v$.
  gsl_vector* SvfromNv0 = gsl_vector_calloc(theta_p);

	// Calculate the forcing term for each time in the period.
  CalcLambda(Lambda, N_v0, eta, theta_p, fntestentopar);

	// Calculate the periodic orbit for the given N_v0.
  CalcXP(x_p, Upsilon, Lambda, inv1Xtp, eta, theta_p, fntestentopar);

	// Extract the number of infectious mosquitoes from the full periodic
	// orbit.
  int indexSv = 2*mt;
  for (size_t i=0; i<theta_p; i++){
    gsl_vector_set(SvfromNv0, i,
                   gsl_vector_get(x_p[i], indexSv));
  }

# ifdef VectorControl_PRINT_CalcSvDiff
  PrintVector(fntestentopar, SvfromNv0name, SvfromNv0, theta_p);
# endif

	// Subtract S_vFromEIR from SvfromNv0
  gsl_vector_memcpy(S_vDiff,SvfromNv0);
  gsl_vector_sub(S_vDiff, S_vFromEIR);


  for (size_t i=0; i<theta_p; i++){
    gsl_vector_free(Lambda[i]);
    gsl_vector_free(x_p[i]);
  }

  free(Lambda);
  free(x_p);
  gsl_vector_free(SvfromNv0);
}


void CalcLambda(gsl_vector** Lambda, gsl_vector* N_v0, int eta,
                int theta_p, char fntestentopar[])
{
  for(size_t t=0; t < theta_p; t++){
    Lambda[t] = gsl_vector_calloc(eta);
    gsl_vector_set(Lambda[t], 0, gsl_vector_get(N_v0, t));
  }
	
# ifdef VectorControl_PRINT_CalcLambda
  // We should try to print some of these vectors out to see what they  look like.
  PrintLambda(Lambda, eta, fntestentopar);		
# endif
}


void CalcXP(gsl_vector** x_p, gsl_matrix** Upsilon, 
            gsl_vector** Lambda, gsl_matrix* inv1Xtp, int eta,
            int theta_p, char fntestentopar[])
{
  gsl_vector* vtemp = gsl_vector_calloc(eta);
	// gsl_vector* vtempsum = gsl_vector_calloc(eta);
  gsl_matrix* mtemp = gsl_matrix_calloc(eta, eta);

	// Initial condition for periodic orbit.
  gsl_vector* x0p = gsl_vector_calloc(eta);

  printf("Entered CalcXP() \n");

	// Evaluate the initial condition of the periodic orbit.
	// Please refer to paper [add reference to paper and equation
	// number here] for the expression for $x_0$.
  for(size_t i=0; i < theta_p; i++){
    FuncX(mtemp, Upsilon, theta_p, i+1, eta);
    gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, Lambda[i], 1.0, vtemp);
		// gsl_vector_add(vtempsum, vtemp);
  }
  gsl_blas_dgemv(CblasNoTrans, 1.0, inv1Xtp, vtemp, 0.0, x0p);

  printf("Calculated initial condition for periodic orbit. \n");
# ifdef VectorControl_PRINT_CalcXP
  char x0pname[15] = "x0p";
  PrintVector(fntestentopar, x0pname, x0p, eta);
# endif

	// We evalute the full periodic orbit now. 
	// Note: to try to keep the indices consistent with our notes and MATLAB, 
	// x_p[0] will refer to x_p(1): because Upsilon[0] refers to Upsilon(1).
	// Thus, x_p[theta_p-1] = x0p. We can check this to make sure.
  for(size_t t=0; t<theta_p; t++){
		// Print t 
    printf("t=%u \r", t+1);
		/*
    if(t==100 || t==200 || t==300){
    printf("t=%d \n", t);
  }
    else{
    printf("t=%d \r", t);
  }
                */
    x_p[t] = gsl_vector_calloc(eta);
		// gsl_vector_set_zero(vtemp);
		// gsl_vector_set_zero(vtempsum);
    FuncX(mtemp, Upsilon, t+1, 0, eta);
    gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, x0p, 1.0, x_p[t]);
    for(size_t i=0; i<=t; i++){
			// printf("t=%d i=%d \n", t, i);
      FuncX(mtemp, Upsilon, t+1, i+1, eta);
      gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, Lambda[i], 1.0, x_p[t]);
    }
  }
        
  printf("Calculated periodic orbit. \n");

# ifdef VectorControl_PRINT_CalcXP
  // We should try to print some of these vectors out to see what they  look like.
  PrintXP(x_p, eta, theta_p, fntestentopar);		
# endif

  gsl_vector_free(vtemp);
	// gsl_vector_free(vtempsum);
  gsl_vector_free(x0p);
  gsl_matrix_free(mtemp);
}


void CalcPSTS(double* sumkplusPtr, double* sumklplus, int theta_s,
              int tau, double P_A, double Pdf)
{
  double taud = (double)tau;
  double thetasd = (double) theta_s;

  int klplus;	// $k_{l+}$ in model.
	// klplus = (int *)malloc((tau-1)*sizeof(int)); Define temporarily.
        // $k_+$ in model:
  int kplus = (int) ((thetasd/taud)-1.); // = floor(theta_s/tau)-1;

	// Evaluate sumkplus
  double sumkplus = 0.;
  for (size_t j=0; j <= kplus; j++){
    double tempbin = binomial(theta_s-(j+1)*tau+j,j);
    double temppap = pow(P_A,(double)theta_s-(j+1)*tau);
    double temppdfp = pow(Pdf,(double)j);
    double temp = tempbin*temppap*temppdfp;
    sumkplus += temp;
  }
  *sumkplusPtr = sumkplus;

	/* // Check results for kplus. 
  itemp = (int) (7./2.) - 1.;
  printf("kplus = %d\n", kplus);
  printf("itemp = %d\n", itemp);
  getchar();
        */

	/* // Print sumkplus
  printf("In CalcPSTS(), sumkplus = %f \n", sumkplus);
  getchar();
        */ 

	// Evaluate sumklplus
  for (size_t l=1; l <= tau-1; l++){
    klplus = (int) (((thetasd+l)/taud) - 2); // = floor((theta_s+l)/tau)-2;
    sumklplus[l-1] = 0;
		// printf("For l = %d, klplus = %d \n", l, klplus);

    for(size_t j=0; j<=klplus; j++){
      double tempbin = binomial(theta_s+l-(j+2)*tau+j,j);
      double temppap = pow(P_A,(double)(theta_s+l-(j+2)*tau));
      double temppdfp = pow(Pdf,(double)(j+1));
      double temp = tempbin*temppap*temppdfp;
      sumklplus[l-1] += temp;
			// printf("For j = %d, tempsum = %f \n", j, temp);
    }
		// printf("sumklplus(%d) = %f \n", l, sumklplus[l-1]);
  }
}


void FuncX(gsl_matrix* X, gsl_matrix** Upsilon, int t, int s, int eta)
{
  gsl_matrix* temp = gsl_matrix_calloc(eta, eta); 

  gsl_matrix_set_identity(X);

  for (size_t i=s; i<t; i++){
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Upsilon[i], X, 0.0, temp);
    gsl_matrix_memcpy(X, temp);
  }
  gsl_matrix_free(temp);
}


double CalcSpectralRadius(gsl_matrix* A, int n, char fntestentopar[])
{
  gsl_vector* abseval = gsl_vector_calloc(n);	// Vector of the absolute values of eigenvalues.
  gsl_matrix* B = gsl_matrix_calloc(n, n); // Use to keep A safe.
  gsl_vector_complex* eval = gsl_vector_complex_calloc(n); // Vector of eigenvalues.
	// Allocate memory for workspace to evaluate the eigenvalues.
  gsl_eigen_nonsymm_workspace* w = gsl_eigen_nonsymm_alloc(n); 

	// Copy A into B to keep it safe.
  gsl_matrix_memcpy(B, A);

	// Calculate eigenvalues of B:
  gsl_eigen_nonsymm(B, eval, w);

# ifdef VectorControl_PRINT_CalcSpectralRadius
  PrintEigenvalues(fntestentopar, eval, n);
# endif

	// Calculate the absolute values of the eigenvalues.
  for(size_t i=0; i<n; i++){
    gsl_complex ztemp = gsl_vector_complex_get(eval, i);
    gsl_vector_set(abseval, i, gsl_complex_abs(ztemp));
  }
	
	// Find the largest eigenvalue.
  double sr;	// sprectral radius
  sr = gsl_vector_max(abseval);

	// Free memory.
  gsl_matrix_free(B);
  gsl_vector_complex_free(eval);
  gsl_eigen_nonsymm_free(w);
  gsl_vector_free(abseval);

  return sr;
}


void CalcInv1minusA(gsl_matrix* inv1A, gsl_matrix* A, int n, char fntestentopar[])
{
	// Data types required to compute inverse.
  gsl_matrix* B = gsl_matrix_calloc(n, n); // We calculate (I-A) in B.
  gsl_permutation* p = gsl_permutation_alloc(n);


  gsl_matrix_set_identity(B); // B = I.
  gsl_matrix_sub(B, A);	// B = I-A.

	// Calculate LU decomposition of (I-A).
  int signum;
  gsl_linalg_LU_decomp(B, p, &signum);

	// Use LU decomposition to calculate inverse.
  gsl_linalg_LU_invert(B, p, inv1A);


# ifdef VectorControl_PRINT_CalcInv1minusA
  char invname[15] = "inv1minusA"; // Name of matrix (when printing to file).
  PrintMatrix(fntestentopar, invname, inv1A, n, n);
# endif

	// Free memory.
  gsl_matrix_free(B);
  gsl_permutation_free(p);
}


void CalSvfromEIRdata(gsl_vector* Sv, double P_Ai, double P_B_i, double N_i, 
                      gsl_vector* Xi_i)
{
  // Sv(t) = Xi_i(t)*(N_i/(P_Ai*P_B_i))
  gsl_vector_memcpy(Sv, Xi_i);
  gsl_vector_scale(Sv, N_i/(P_Ai*P_B_i));	
}


double binomial(int n, int k)
{
  return gsl_sf_fact((unsigned int) n) / (gsl_sf_fact((unsigned int) k) * gsl_sf_fact((unsigned int) (n-k)));
}



/******************************************************************************
  Printing routines below. Most are only optionally compiled in.
******************************************************************************/

void PrintRootFindingStateTS(size_t iter, gsl_multiroot_fsolver* srootfind, 
                             int theta_p, char fnrootfindingstate[])
{
/* PrintRootFindingStateTS() prints the current status of the root-
  * finding algorithm to the screen and to the given file.
  *
  * There are numerous quantities that we could print to see how the
  * root-finding algorithm is doing. It is not reasonable to print
  * all theta_p terms, so for now, we print out the value of N_v0[0]
  * to see one of the values of the emergence rate, and the $l^1$
  * norm of $f$.
  *
  * Note that we print to screen and to a file.
  *
  * All parameters are IN parameters.
 */   
  double svdiffsum;
  double Nv0_0;

  FILE* fpp = fopen(fnrootfindingstate, "a");

	// Calculate the $l^1$ norm of f.
  svdiffsum = gsl_blas_dasum(srootfind->f);

	// Get the 0th element of N_v0.
  Nv0_0 = gsl_vector_get(srootfind->x, 0);

	// Print to screen:
  printf("iter = %5u N_v0(1) = % .3f ||f||_1 = % .3f \n", iter, Nv0_0, svdiffsum);
  fprintf(fpp, "iter = %5u N_v0(1) = % .3f ||f||_1 = % .3f \n", iter, Nv0_0, svdiffsum);
  fclose(fpp);
}

#ifdef VectorControl_PRINT_CalcInitMosqEmergeRate	// only use
void PrintParameters(char fntestentopar[], int theta_p, int tau, int theta_s, 
                    int n, int m, double N_i, double alpha_i, double mu_vA, 
                    double theta_d, double P_B_i, double P_C_i, double P_D_i, double P_E_i, 
                    gsl_vector* K_vi, gsl_vector* Xi_i, gsl_vector* Nv0guess)
{
/* PrintParameters() prints the input parameters to a given file. 
  * We currently use this to make sure that the inputs we have in C
  * are what we expect from what we've sent from Fortran. 
  *
  * We may transform/copy this into a new function that does more.
  * 
  * All parameters are IN parameters.
 */
  int i;
  double temp;

  FILE* fpp = fopen(fntestentopar, "a");

  fprintf(fpp, "theta_p = %d; \n", theta_p);
  fprintf(fpp, "tau = %d; \n", tau);
  fprintf(fpp, "theta_s = %d; \n", theta_s);
  fprintf(fpp, "n = %d; \n", n);
  fprintf(fpp, "m = %d; \n", m);

  fprintf(fpp, "N_i = %f; \n", N_i);
  fprintf(fpp, "alpha_i = %f; \n", alpha_i);
  fprintf(fpp, "mu_vA = %f; \n", mu_vA);
  fprintf(fpp, "theta_d = %f; \n", theta_d);
  fprintf(fpp, "P_B_i = %f; \n", P_B_i);
  fprintf(fpp, "P_C_i = %f; \n", P_C_i);
  fprintf(fpp, "P_D_i = %f; \n", P_D_i);
  fprintf(fpp, "P_E_i = %f; \n", P_E_i);

	
  fprintf(fpp, "K_vi = \n");
  gsl_vector_fprintf(fpp, K_vi, "%f");

  fprintf(fpp, "Xi_i = \n");
  gsl_vector_fprintf(fpp, Xi_i, "%f");

  fprintf(fpp, "Nv0guess = \n");
  gsl_vector_fprintf(fpp, Nv0guess, "%f");
	
	// Let's do this properly.
	
  for (i=0; i<theta_p; i++){
    temp = gsl_vector_get(K_vi, i);
    fprintf(fpp, "K_vi(%d) = %f; \n", i+1, temp);
  }

  for (i=0; i<theta_p; i++){
    temp = gsl_vector_get(Xi_i,i);
    fprintf(fpp, "Xi_i(%d) = %f; \n", i+1, temp);
  }

  for (i=0; i<theta_p; i++){
    temp = gsl_vector_get(Nv0guess, i);
    fprintf(fpp, "Nv0guess(%d) = %f; \n", i+1, temp);
  }
	
  fclose(fpp);
}
#endif

#ifdef VectorControl_PRINT_CalcUpsilonOneHost	// only use
void PrintUpsilon(char fntestentopar[], gsl_matrix** Upsilon, int theta_p,
                  int eta, double P_A, double P_Ai, double Pdf, gsl_vector* Pdif,
                  gsl_vector* Pduf)
{
/* PrintUpsilon() prints the intermediate results while calculating 
  * Upsilon.
  * 
  * All parameters are IN parameters.
 */

  int i;
  int j;
  int k;
  double temp;

  FILE* fpp = fopen(fntestentopar, "a");

  fprintf(fpp, "P_A = %f\n", P_A);
  fprintf(fpp, "P_Ai = %f\n", P_Ai);
  fprintf(fpp, "Pdf = %f\n", Pdf);

	/*
  for (i=0; i<theta_p; i++){
  temp = gsl_vector_get(Pdif, i);
  fprintf(fpp, "Pdif(%d) = %f \n", i+1, temp);
}

  for (i=0; i<theta_p; i++){
  temp = gsl_vector_get(Pduf, i);
  fprintf(fpp, "Pduf(%d) = %f \n", i+1, temp);
}
        */ 

	// Print some Upsilon[k].
  k=0;

  fprintf(fpp, "Upsilon[%d] = \n", k);
  for (i=0; i < eta; i++){
    for (j=0; j < eta; j++){
      temp = gsl_matrix_get(Upsilon[k], i, j);
      fprintf(fpp, "%f ", temp);
    }
    fprintf(fpp, "\n");
  }
	// gsl_matrix_fprintf(fpp, Upsilon[k], "%f");

  k = 364;

  fprintf(fpp, "Upsilon[%d] = \n", k);
  for (i=0; i < eta; i++){
    for (j=0; j < eta; j++){
      temp = gsl_matrix_get(Upsilon[k], i, j);
      fprintf(fpp, "%f ", temp);
    }
    fprintf(fpp, "\n");
  }

  fclose(fpp);
}
#endif

#ifdef VectorControl_PRINT_CalcXP	// only use
                                                                            void PrintXP(gsl_vector** x_p, int eta, int theta_p, char fntestentopar[])
{
/* PrintXP() prints out values of XP, the periodic orbit.
  * 
  * All parameters are IN parameters.
 */
  int t;
  char xpname[15] = "x_p";
  char xpvecname[15];
  char timestring[5];
	// double temp;

  FILE* fpp = fopen(fntestentopar, "a");

	// Print all x_p[t]:
  for(t=0; t<theta_p; t++){
    sprintf(timestring, "%d", t+1);
    strcpy(xpvecname, xpname);
    strcat(xpvecname, "(");
    strcat(xpvecname, timestring);
    strcat(xpvecname, ")");
    PrintVector(fntestentopar, xpvecname, x_p[t], eta);
  }

  fclose(fpp);
}
#endif

#ifdef VectorControl_PRINT_CalcLambda	// only use
                                                                            void PrintLambda(gsl_vector** Lambda, int eta, char fntestentopar[])
{
/* PrintLambda() prints some values of Lambda.
  * 
  * All parameters are IN parameters.
 */
  int t;
	// double temp;

  FILE* fpp = fopen(fntestentopar, "a");

	// Print some Lambda[t].
  t=0;

  fprintf(fpp, "Lambda[%d] = \n", t);
  gsl_vector_fprintf(fpp, Lambda[t], "%f");

  t=139;

  fprintf(fpp, "Lambda[%d] = \n", t);
  gsl_vector_fprintf(fpp, Lambda[t], "%f");

  t=363;

  fprintf(fpp, "Lambda[%d] = \n", t);
  gsl_vector_fprintf(fpp, Lambda[t], "%f");

  fclose(fpp);
}
#endif

#ifdef VectorControl_PRINT_CalcSpectralRadius	// only use
                                                                            void PrintEigenvalues(char fntestentopar[],gsl_vector_complex* eval, int n)
{
/* PrintEigenvalues() prints eigenvalues to the given file.
  * 
  * All parameters are IN parameters.
 */
	// int i;
	// int j;
	// double temp;

  FILE* fpp = fopen(fntestentopar, "a");

  fprintf(fpp, "Eigenvalues = \n");
	/*
  for (i=0; i < eta; i++){
  for (j=0; j < eta; j++){
  temp = gsl_matrix_get(X_t_p, i, j);
  fprintf(fpp, "%e ", temp);
}
  fprintf(fpp, "\n");
}
        */
  gsl_vector_complex_fprintf(fpp, eval, "%e");
  fclose(fpp);
}
#endif

#if defined VectorControl_PRINT_CalcInitMosqEmergeRate || defined VectorControl_PRINT_CalcInv1minusA
void PrintMatrix(char fntestentopar[], char matrixname[], gsl_matrix* A, 
                 int RowLength, int ColLength)
{

/* PrintMatrix() prints the given matrix to the given file.
  * 
  * All parameters are IN parameters.
 */

  int i;
  int j;
  double temp;

  FILE* fpp = fopen(fntestentopar, "a");

  fprintf(fpp, "%s = \n", matrixname);
  for (i=0; i < ColLength; i++){
    for (j=0; j < RowLength; j++){
      temp = gsl_matrix_get(A, i, j);
      fprintf(fpp, "%e ", temp);
    }
    fprintf(fpp, "\n");
  }

  fclose(fpp);
}
#endif

#if defined VectorControl_PRINT_CalcInitMosqEmergeRate || defined VectorControl_PRINT_CalcSvDiff || defined VectorControl_PRINT_CalcXP
void PrintVector(char fntestentopar[], char vectorname[], gsl_vector* v, int n)
{
/* PrintVector() prints the given (GSL) vector to the given file.
  * 
  * All parameters are IN parameters.
 */
  int i;
  double temp;

  FILE* fpp = fopen(fntestentopar, "a");

  for (i=0; i < n; i++){
    temp = gsl_vector_get(v, i);
    fprintf(fpp, "%s(%d) = %f; \n", vectorname, i+1, temp);
		
  }
	
  fclose(fpp);
}
#endif
