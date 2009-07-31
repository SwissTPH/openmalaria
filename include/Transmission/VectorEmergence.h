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

/// This is a header for VectorTransmission's functions which are only used internally.

#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf.h>
#include <fstream>
#include <limits>
#include "global.h"
using namespace std;

#define VectorTransmission_PRINT_CalcInitMosqEmergeRate
#define VectorTransmission_PRINT_CalcUpsilonOneHost
#define VectorTransmission_PRINT_CalcSvDiff
// #define VectorTransmission_PRINT_CalcLambda
// #define VectorTransmission_PRINT_CalcXP
#define VectorTransmission_PRINT_CalcSpectralRadius
#define VectorTransmission_PRINT_CalcInv1minusA


/** Container for functions used to calculate the mosquito emergence rate.
 *
 * Some data is stored here, so that it doesn't have to be continually
 * freed and reallocated. It is cleaned by the destructor.
 * 
 * All non-const data outside of functions should be stored in here, so as
 * to be thread-safe.
 * 
 * Const-ness: Many of the params passed in to functions are now const.
 * Making functions themselves const seems utterly pointless, since in C++ in
 * const functions, only values directly stored in the object are effectively
 * const, and anything pointed by them (including array members) is non-const.
 */
class VectorEmergence {
public:
  /** Initialises some data elements */
  VectorEmergence(int mosqRestDuration, int EIPDuration, int populationSize,
		  double entoAvailability, double mosqSeekingDeathRate, double mosqSeekingDuration,
		  double probMosqBiting, double probMosqFindRestSite, double probMosqSurvivalResting,
		  double probMosqSurvivalOvipositing,
		  int yearLength = daysInYear, ostream& traceOut = cout,
		  const char* logFileName = "output_ento_para.txt");
  /** Frees data */
  ~VectorEmergence();
  
  /** calcInitMosqEmergeRate() calculates the mosquito emergence rate given
   * all other parameters.
   *
   * We use a periodic version of the model described in "A Mathematical Model 
   * for the Dynamics of Malaria in Mosquitoes Feeding on a Heteregeneous Host
   * Population". The periodic model still needs to be written as a paper. We will
   * change these comments to refer to the approprirate paper when it is ready.
   *
   * The entomological model has a number of input parameters, including the
   * mosquito emergence rate, \f$N_{v0}\f$, and a number of output parameters, 
   * including the entomological inoculation rate, \f$\Xi_i\f$. The model produces
   * equations for \f$\Xi_i\f$ as a function of \f$N_{v0}\f$ and the other parameters.
   * However, in this function, we assume that all parameters, except \f$N_{v0}\f$ 
   * are known, and \f$\Xi_i\f$ is known. We then use these parameters, with \f$\Xi_i\f$ 
   * to calculate \f$N_{v0}\f$. The equations for \f$\Xi_i\f$ are linear in terms of 
   * \f$N_{v0}\f$ so there is a unique solution for \f$N_{v0}\f$. 
   *
   * This routine first shows the existence of a unique globally asymptotically 
   * stable periodic orbit for the system of equations describing the periodically
   * forced entomological model (for a given set of parameter values, including the
   * mosquito emergence rate). It then compares the number of infectious host-seeking
   * mosquitoes for this periodic orbit to the the number of infectious host-seeking
   * mosquitoes that would result in the given EIR. The routine then iteratively finds
   * the emergence rate that matches the given EIR.
   * 
   * However, we cannot write these equations in the form Ax=b, so we use
   * a root-finding algorithm to calculate \f$N_{v0}\f$.
   *
   * This function has a dummy return of 0.
   * 
   * \param mosqEmergeRate is both an input (guessed or read from file) and output (calculated emergence rate-. */
  double CalcInitMosqEmergeRate(int nHostTypesInit, int nMalHostTypesInit,
                                const double* FHumanInfectivityInitVector,
                                const vector<double>& FEIRInitVector,
				double* mosqEmergeRate);
  
private:
  //BEGIN data
  // n and m from the model are not renamed such here; they are:
  // nHostTypesInit, nMalHostTypesInit
  
  ///@brief Parameters that help to describe the order of the system.
  //@{
  /// Ask not why we call mt, mt. We use mt to index the system.
  /// It is the maximum number of time steps we go back for \f$N_v\f$ and \f$O_v\f$.
  size_t mt;
  /// \f$\eta\f$: The order of the system.
  size_t eta;
  //@}
  
  size_t counterSvDiff;
  size_t theta_p;
  size_t tau;
  size_t theta_s;
  
  int N_i;
  double alpha_i;
  double mu_vA;
  double theta_d;
  double P_B_i;
  double P_C_i;
  double P_D_i;
  double P_E_i;
  
  /** The set of theta_p matrices that determine the dynamics of the system
   * from one step to the next.
   *
   * That is, the system is described by,
   * \f$x(t) = \Upsilon(t) x(t-1) = \Lambda(t)\f$.
   * \f$\Upsilon(t)\f$ is defined over time, \f$1 \leq t \leq \theta_p\f$, 
   * where \f$t \in \mathbb{N}\f$. */
  gsl_matrix** Upsilon;
  
  /** The set of theta_p vectors that determine the forcing of the system
   * at every time step.
   *
   * \f$\Lambda(t)\f$ is defined over time, \f$1 \leq t \leq \theta_p\f$, 
   * where \f$t \in \mathbb{N}\f$. */
  gsl_vector** Lambda;
  
  /// The periodic orbit of all eta state variables.
  gsl_vector** x_p;
  
  /// @brief Cached memory; values only have meaning within some functions.
  //@{
  gsl_vector* memVectorEta;	///< eta long vector
  gsl_vector_complex* memVectorComplexEta;	///< eta long vector of complex values
  gsl_matrix* memMatrixEtaSq;	///< eta by eta matrix
  gsl_eigen_nonsymm_workspace* memEigenWorkspace;	///< for dimension eta
  gsl_permutation* memPermutation;	///< dimension eta
  //@}
  
  ostream& trace;	///< Run-time output is printed here
  mutable ofstream logFile;	///< Values are logged to here
  //END data
  
/** CalcUpsilonOneHost returns a pointer to an array of theta_p 
 * GSL matrices assuming there is only one host of humans..
 * Each matrix is Upsilon(t).
 *
 * \f$Upsilon(t)\f$ is the evolution of the mosquito population over one
 * time step. There are three main system variables:
 * \f$N_v\f$: The total number of host-seeking mosquitoes.
 * \f$O_v\f$: The number of infected host-seeking mosquitoes.
 * \f$S_v\f$: The number of infectious host-seeking mosquitoes.
 *
 * As the difference equations go back more than one time step, 
 * the size of the system is larger than 3.
 * For \f$N_v\f$ and \f$O_v\f$, we need to go back mt steps.
 * For \f$S_v\f$ we need to go back tau steps.
 * So the size of the system, eta = 2 mt + tau.
 * The first column of Upsilon(t) (indexed by 0 in C) corresponds to
 * \f$N_v(t)\f$ - as it depends on the other paramters at previous times.
 * The (mt+1)^th column of Upsilon(t) (indexed by mt in C) corresponds to
 * \f$O_v(t)\f$ - as it depends on the other paramters at previous times.
 * The (2mt+1)^th column of Upsilon(t) (indexed by 2mt in C) corresponds to
 * \f$S_v(t)\f$ - as it depends on the other paramters at previous times.
 * All other columns have 1 in the subdiagonal.
 *
 * For now, we write this code assuming that the parameters where we
 * are ignoring dependence on host type, or phase of the period, (and
 * have defined as doubles) will remove doubles. We do not code for
 * generality. If we make changes to these data types later, we will
 * change the code then. We code this, as is, to make it easier now, 
 * as we do not know what parameters we will change. It should 
 * hopefully, not be too difficult to change the code later (and
 * create a new general CalcUpsilon). Let's hope....
 *
 * Upsilon is set.
 * 
 * PAPtr, and PAiPtr are OUT parameters.
 * All other parameters are IN parameters. */
void CalcUpsilonOneHost(double* PAPtr, double* PAiPtr, size_t n, size_t m, const gsl_vector* K_vi);

/** CalcSvDiff returns the difference between Sv for the periodic 
 * orbit for the given Nv0 and from the EIR data.
 * 
 * Given the input parameters to the entomological model, this routine
 * calculates the number of infectious host-seekign mosquitoes for the 
 * resulting periodic orbit. It then calculates the difference between 
 * this Sv and the periodic Sv calculated from the EIR data (which is 
 * the Sv from the periodic orbit of the system with the final 
 * calculated Nv0.
 * 
 * Upsilon is read.
 * 
 * SvDiff is an OUT parameter.
 * All other parameters are IN parameters. */
void CalcSvDiff(gsl_vector* SvDiff, const gsl_vector* SvfromEIR,
                const gsl_vector* Nv0, const gsl_matrix* inv1Xtp);

/** CalcLambda() returns a pointer to an array of theta_p 
 * GSL vectors.
 * Each vector is Lambda(t).
 *
 * \f$Lambda(t)\f$ is the forcing of the mosquito population
 * at each time step, that is, it is the number of new
 * mosquitoes that enter the population at each time, \f$t\f$.
 *
 * We note here that Nv0 is a gsl_vector where the index, t, refers
 * to the mosquito emergence rate at time, t. Lambda[t] is a
 * gsl_vector that denotes the forcing at time t, where the index, i,
 * refers to the forcing to the i^th dimension of the system.
 *
 * God moves over the face of the waters,
 * Looking to the left and looking to the right,
 * But there is only water to see.
 * 
 * Lambda is set.
 * 
 * All parameters are IN parameters. */
void CalcLambda(const gsl_vector* Nv0);

/** CalcXP returns a pointer to an array of theta_p 
 * GSL vectors.
 * Each vector is is the periodic orbit solution to the main system
 * of equations at time, t.
 *
 * The size of each xp[t] is eta: the order of the system.
 *
 * This routine uses Theorem 2 of Cushing (1998) JDEA 3.
 *
 * We can probably improve the speed of this algorithm. We could replace
 * the vectors for Lambda[i] by simply using Nv0[i] and multiplying the first
 * column of the matrices, X[t,i] by Nv0[i]. 
 *
 * But for now we don't worry about speed and try to continue with the route
 * finding. There may be more work that we need to do to improve speed.
 * 
 * Upsilon, Lambda are read.
 * x_p is set.
 * 
 * All parameters are IN parameters. */
void CalcXP(const gsl_matrix* inv1Xtp);


/** CalcPSTS() calculates probabilities of surviving the extrinsic
 * incubation period (or part of). The returned variables are the sums
 * to \f$k_+\f$ and \f$k_{l+}\f$ (including the binomial coefficients and 
 * probabilities in (2.3c) of the paper. 
 *
 * Currently, this returns scalar values because neither \f$P_A\f$, nor
 * \f$P_{df}\f$, depend on the phase of the period.
 *
 * Note that sumklplus here is defined as sumlv in MATLAB.
 * 
 * sumkplusPtr and sumklplus are OUT parameters.
 * All other parameters are IN parameter. */
void CalcPSTS(double* sumkplusPtr, double* sumklplus, double P_A, double P_df) const;

/** FuncX() calculates X(t,s).
 *
 * Note that we have to be careful with indices here. 
 * Cushing (1995) has indices starting at 0 and ending at \f$\theta_p -1\f$.
 * In our notes, and in MATLAB, the indices start at 1 and end at \f$\theta_p\f$.
 *
 *       X(t,s) = \Upsilon(t-1)*...*Upsilon(s) for t \geq s+1
 *              = I                            for t = s.
 *
 * Here, FuncX() is defined for s>=0 and t>=1.
 * 
 * Upsilon is read.
 * 
 * X is an OUT parameter, t and s are IN parameters. */
void FuncX(gsl_matrix* X, size_t t, size_t s) const;

/** CalcSpectralRadius() calculates the spectral radius of a given matrix.
 *
 * Given an eta by eta, real, nonsymmetric matrix, A, 
 * this routine calcultes its spectral radius,
 * that is, the eigenvalue with the largest absolute value.
 * 
 * Dimensions are required to be eta×eta in order to allow re-using memory.
 * 
 * A is an IN parameters. */
double CalcSpectralRadius(const gsl_matrix* A) const;

/** CalcInv1minusA() calculates the inverse of (I-A) where A is a 
 * given matrix.
 *
 * Given an eta by eta, real matrix, A, 
 * this routine calcultes the inverse of (I-A) where I is the 
 * eta by eta identity matrix.
 * 
 * Dimensions are required to be eta×eta in order to allow re-using memory.
 * 
 * A is an IN parameter.
 * inv1A is an OUT parameter. */
void CalcInv1minusA(gsl_matrix* inv1A, const gsl_matrix* A) const;

/** CalcSvfromEIRdata() calculates Sv, given the EIR.
 *
 * Given EIR, and the parameters that determine host-biting,
 * this routine calculates the number of infectious host-seeking
 * mosquitoes, Sv.
 *
 * The EIR is assumed to be periodic so the resulting vector for Sv
 * is also periodic.
 *
 * The other parameters are constant. 
 *
 * Once the periodic paper is written, we should add a refernece to the 
 * equation that we use.
 * 
 * P_Ai, P_B_i and Xi_i are IN parameters.
 * Sv is an OUT parameter. */
void CalSvfromEIRdata(gsl_vector* Sv, double P_Ai, const gsl_vector* Xi_i) const;


/** binomial() calculates the binomial coefficient given two integers.
 * 
 * Note that we do not check for errors. */
double binomial(int n, int k) const;


/******************************************************************************
 Printing routines below. Most are only optionally compiled in.
******************************************************************************/
/** PrintRootFindingStateTS() prints the current status of the root-
  * finding algorithm to the screen and to the given file.
  *
  * There are numerous quantities that we could print to see how the
  * root-finding algorithm is doing. It is not reasonable to print
  * all theta_p terms, so for now, we print out the value of N_v0[0]
  * to see one of the values of the emergence rate, and the \f$l^1\f$
  * norm of \f$f\f$.
  *
  * Note that we print to screen and to a file.
  *
  * All parameters are IN parameters.
 */
void PrintRootFindingStateTS(size_t iter, const gsl_multiroot_fsolver* srootfind) const;

/** PrintParameters() prints the input parameters to a given file. 
  * We currently use this to make sure that the inputs we have in C
  * are what we expect from what we've sent from Fortran. 
  *
  * We may transform/copy this into a new function that does more.
  * 
  * All parameters are IN parameters.
 */
void PrintParameters(size_t theta_p, size_t tau, size_t theta_s,
		size_t n, size_t m, double N_i, double alpha_i, double mu_vA,
		double theta_d, double P_B_i, double P_C_i, double P_D_i, double P_E_i,
		const gsl_vector* K_vi, const gsl_vector* Xi_i) const;

/** PrintUpsilon() prints the intermediate results while calculating 
  * Upsilon.
  * 
  * All parameters are IN parameters.
 */
void PrintUpsilon(const gsl_matrix *const * Upsilon, size_t theta_p,
		size_t eta, double P_A, double P_Ai, double P_df,
		const gsl_vector* P_dif, const gsl_vector* P_duf) const;

/** PrintXP() prints out values of XP, the periodic orbit.
  * 
  * All parameters are IN parameters.
 */
void PrintXP(const gsl_vector *const * x_p, size_t eta, size_t theta_p) const;

void PrintMatrix(const char matrixname[], const gsl_matrix* A,
		size_t RowLength, size_t ColLength) const;

public:
  /** PrintVector() prints the given (GSL) vector to the given file. */
  void PrintVector(const char* vectorname, const gsl_vector* v) const;

  /** PrintArray() prints the given (C) array to the given file.
   * 
   * The array, v, of doubles is assumed to be of length n. */
  void PrintArray(const char* vectorname, const double* v, int n) const;
  /// ditto, taking a vector
  void PrintArray(const char* vectorname, const vector<double>& v) const;

  friend int CalcSvDiff_rf(const gsl_vector* x, void* p, gsl_vector* f);
  friend class VectorEmergenceSuite;
};

/** @brief Free functions called by the GSL root finder
 *
 * GSL can't take member-function pointers, so instead we pass free functions,
 * with a pointer to the VectorEmergence object. */
//@{
/** CalcSvDiff_rf returns the difference between Sv for the periodic 
 * orbit for the given Nv0 and from the EIR data.
 * 
 * Given the input parameters to the entomological model, this routine
 * calculates the number of infectious host-seekign mosquitoes for the 
 * resulting periodic orbit. It then calculates the difference between 
 * this Sv and the periodic Sv calculated from the EIR data (which is 
 * the Sv from the periodic orbit of the system with the final 
 * calculated Nv0.
 *
 * This routine performs the same calculations as CalcSvDiff() but it
 * does so in the format required by the GSL multiroot-finding
 * algorithms.
 * 
 * f is an OUT parameter.
 * All other parameters are IN parameters. */
int CalcSvDiff_rf(const gsl_vector* x, void* p, gsl_vector* f);
//@}

// Structure that contains the parameters for the function used in the 
// root-finding algorithm to find the emergence rate that matches the 
// number of infectious host-seeking mosquitoes.
struct SvDiffParams
{	//FIXME: move some of these to VectorEmergence
  SvDiffParams (VectorEmergence* e, gsl_vector* v, gsl_matrix* m, size_t theta_p) :
    emerge(e), S_vFromEIR(v), inv1Xtp(m),
    lastNv0(gsl_vector_calloc (theta_p)), lastS_vDiff(gsl_vector_calloc (theta_p))
  {
    // make sure lastNv0 won't match any input first time, so lastS_vDiff will be calculated
    gsl_vector_set (lastNv0, 0, numeric_limits<double>::quiet_NaN());
  }
  ~SvDiffParams () {
    gsl_vector_free (lastNv0);
    gsl_vector_free (lastS_vDiff);
  }
  
  VectorEmergence* emerge;
  gsl_vector* S_vFromEIR;
  gsl_matrix* inv1Xtp;
  /// The last Nv0 vector used to calculate x_p and S_vDiff
  /// (if not current, we recalculate x_p during the root finding).
  gsl_vector* lastNv0;
  /// The last calculated S_vDiff
  gsl_vector* lastS_vDiff;
};
