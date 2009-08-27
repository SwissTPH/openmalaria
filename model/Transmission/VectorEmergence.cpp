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


/* --- This file contains free functions used by VectorTransmission  ---
   ---         (to calculate mosquito emergence rates).         --- */


#include "Transmission/VectorEmergence.h"
#include "global.h"
#include "util/vectors.h"

#include <sstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <string.h>

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
VectorEmergence::VectorEmergence(int mosqRestDuration, int EIPDuration,
    int populationSize, double entoAvailability, double mosqSeekingDeathRate,
    double mosqSeekingDuration, double probMosqBiting, double probMosqFindRestSite,
    double probMosqSurvivalResting, double probMosqSurvivalOvipositing,
    int yearLength, ostream& traceOut, const char* logFileName) :
  counterSvDiff(0),
  theta_p (yearLength), tau(mosqRestDuration), theta_s(EIPDuration),
  N_i(populationSize),
  alpha_i(entoAvailability),
  mu_vA(mosqSeekingDeathRate), theta_d(mosqSeekingDuration),
  P_B_i(probMosqBiting), P_C_i(probMosqFindRestSite),
  P_D_i(probMosqSurvivalResting), P_E_i(probMosqSurvivalOvipositing),
  trace (traceOut), logFile (logFileName)
{
  mt = theta_s + tau -1;
  eta = 2*mt + tau;
  
  Upsilon = new gsl_matrix*[theta_p];
  Lambda = new gsl_vector*[theta_p];
  x_p = new gsl_vector*[theta_p];
  for (size_t i=0; i<theta_p; ++i) {
    Upsilon[i] = gsl_matrix_calloc(eta, eta);
    Lambda[i] = gsl_vector_calloc(eta);
    x_p[i] = gsl_vector_calloc(eta);
  }
  
  memVectorEta = gsl_vector_calloc(eta);
  memVectorComplexEta = gsl_vector_complex_calloc(eta);
  memMatrixEtaSq = gsl_matrix_calloc(eta, eta);
  memEigenWorkspace = gsl_eigen_nonsymm_alloc(eta);
  memPermutation = gsl_permutation_alloc(eta);
}

VectorEmergence::~VectorEmergence () {
  gsl_vector_complex_free(memVectorComplexEta);
  gsl_vector_free(memVectorEta);
  gsl_matrix_free(memMatrixEtaSq);
  gsl_eigen_nonsymm_free(memEigenWorkspace);
  gsl_permutation_free(memPermutation);
  
  for (size_t i=0; i<theta_p; i++){
    gsl_matrix_free(Upsilon[i]);
    gsl_vector_free(Lambda[i]);
    gsl_vector_free(x_p[i]);
  }
  delete[] Upsilon;
  delete[] Lambda;
  delete[] x_p;
  
  logFile.close();
}

bool VectorEmergence::CalcInitMosqEmergeRate(const vector<double>& FHumanInfectivityInitVector,
					     const vector<double>& FEIRInitVector,
					     const NonHumanHostsType& nonHumanHosts,
					     vector<double>& mosqEmergeRate)
{
  PrintVector ("kappa", FHumanInfectivityInitVector);
  
  // This initially contains the initial estimate of the mosquito emergence
  // rate. This is used by the root finding algorithm to calculate N_v0.
  gsl_vector* N_v0 = vectors::std2gsl (mosqEmergeRate, theta_p);
  PrintVector("Nv0", N_v0);
  
  gsl_vector* K_vi = vectors::std2gsl (FHumanInfectivityInitVector, theta_p);
  
  // Output Parameters (for the model):
  gsl_vector* Xi_i = vectors::std2gsl (FEIRInitVector, theta_p);
  
  // The number of infectious mosquitoes over every day of the cycle.
  // calculated from the EIR data.
  // \f$S_v\f$ (from EIR).
  gsl_vector* S_vFromEIR = gsl_vector_calloc(theta_p);
  
  // The difference between S_vFromEIR and SvfromNv0.
  gsl_vector* S_vDiff = gsl_vector_calloc(theta_p);
  
  // State variables.
  
  // Allocate memory for gsl_matrices and initialize to 0.
  
  // \f$X_{\theta_p}\f$.
  // The product of all the evolution matrices.
  // \f$X_{\theta_p} = X(\theta_p+1,1)\f$. 
  // Refer to Cushing (1995) and the paper for the periodic entomological model
  // for more information.
  gsl_matrix* X_t_p = gsl_matrix_calloc(eta, eta);
  
  // \f$(\mathbb{I}-X_{\theta_p})^{-1}\f$.
  // The inverse of the identity matrix minus X_t_p.
  gsl_matrix* inv1Xtp = gsl_matrix_calloc(eta, eta);

# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  // We now try to print these parameters to file to make sure that 
  // they show what we want them to show.
  PrintParameters(theta_p, tau, theta_s, N_i, alpha_i,
                  mu_vA, theta_d, P_B_i, P_C_i, P_D_i, P_E_i, K_vi, Xi_i);
  // The parameter values look correct.
# endif
  
  // Derived Parameters
  
  /// Probability that a mosquito survives one day of 
  /// host-seeking but does not find a host. 
  double P_A = 0;

  /// Probability that on a given day, a mosquito finds a host of type \f$i\f$.
  /// For now, we assume that this  is  a double: 
  /// - no dependence on the phase of the period - or the type of host.
  double P_Ai = 0;	// \f$P_{A^i}\f$
  
  // Create matrices in Upsilon.
  // We also define P_A and P_Ai in the same routine. 
  // For now, we treat P_A and P_Ai as scalars since we are 
  // defining most parameters as scalars. If we do change things later, which we
  // may, then we will change the code accordingly. We will need to go through
  // a lot of changes anyway. 
  CalcUpsilonOneHost(&P_A, &P_Ai, K_vi, nonHumanHosts);
  
  // Calculate \f$X_{\theta_p}\f$.
  // Refer to Cushing (1995) and the paper for the periodic entomological model
  // for more information.
  FuncX(X_t_p, theta_p, 0);
  
# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  char xtpname[15] = "X_t_p";
  PrintMatrix(xtpname, X_t_p, eta, eta);
# endif

  // We should now find the spectral radius of X_t_p and show that it's less than 1.
  double srXtp = CalcSpectralRadius(X_t_p);

  // trace << "The spectral radius of X_t_p = " << srXtp << endl;
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
    msg << "Warning: All results from the entomologoical model may be meaningless.";
    throw xml_scenario_error (msg.str());
  }
  
  // Calculate the inverse of (I-X_t_p).
  CalcInv1minusA(inv1Xtp, X_t_p);

# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  char inv1Xtpname[15] = "inv1minusXtp";
  PrintMatrix(inv1Xtpname, inv1Xtp, eta, eta);
# endif
  
  // Calculate the number of infectious host-seeking mosquitoes for the given EIR.
  CalSvfromEIRdata(S_vFromEIR, P_Ai, Xi_i);

# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  PrintVector("S_vFromEIR", S_vFromEIR);
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
  
  
  CalcSvDiff(S_vDiff, S_vFromEIR, N_v0, inv1Xtp);
# ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
  PrintVector("InitSvDifference", S_vDiff);
# endif
  double SvDiff1norm = gsl_blas_dasum(S_vDiff);	//The $l^1$ norm of S_vDiff.
  trace << "The $l^1$ norm of S_vDiff is " << SvDiff1norm << endl;

  // Maximum \f$l^1\f$ distance of error of root-finding algorithm
  const double EpsAbsRF = 1.0;
  
  bool valid = SvDiff1norm <= EpsAbsRF;
  if (!valid){
    trace << "The difference in Sv is greater than the tolerance." << endl;
    
    if (!(Global::clOptions & CLO::ENABLE_ERC))
      throw runtime_error ("Cannot recalculate: emergence rate calculations are not enabled.");
    
    /************* We initialize variables for root-finding. **************/
    trace << "Starting root-finding" << endl;

    // Parameters for root-finding function.
    struct SvDiffParams pararootfind (this, S_vFromEIR, inv1Xtp,
				      theta_p);
    
    // Set root-finding function.
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

    trace << "About to set root-finding solver" << endl;
    // Initialize root-finding.
    gsl_multiroot_fsolver_set(srootfind, &frootfind, xrootfind);
    trace << "Set root-finding" << endl;

    // Print initial state (to screen):
    PrintRootFindingStateTS(0, srootfind);
    
    // Maximum number of iterations of root-finding algorithm.
    const size_t maxiterRF = 1000;
    int status = GSL_CONTINUE;
    for (size_t iter = 1; status == GSL_CONTINUE && iter < maxiterRF; ++iter) {
      status = gsl_multiroot_fsolver_iterate(srootfind);
      PrintRootFindingStateTS(iter, srootfind);

      // Check to see if solver is stuck
      if (status){
        break;
      }

      status = gsl_multiroot_test_residual(srootfind->f, EpsAbsRF);
    }

    // Print status
    trace << "status = " <<  gsl_strerror(status) << endl; 

    // Copy solution for N_v0 into N_v0.
    gsl_vector_memcpy(N_v0, srootfind->x);

#   ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate
    PrintVector("FinalNv0", N_v0);
    PrintVector("FinalSvDiff", srootfind->f);
#   endif

		// Free memory.
    gsl_vector_free(xrootfind);
    gsl_multiroot_fsolver_free(srootfind);
  }
  

  // Calculate final periodic orbit.
  CalcLambda(N_v0);
  CalcXP(inv1Xtp);
  
  // Copy the mosquito emergence rate to the C array.
  memcpy (&mosqEmergeRate[0], N_v0->data, theta_p * sizeof (double));
  
  gsl_vector_free(N_v0);
  gsl_vector_free(K_vi);
  gsl_vector_free(Xi_i);
  gsl_vector_free(S_vFromEIR);
  gsl_vector_free(S_vDiff);
  
  gsl_matrix_free(X_t_p);
  gsl_matrix_free(inv1Xtp);
  
  return valid;
}


void VectorEmergence::CalcUpsilonOneHost(double* PAPtr, double* PAiPtr, const gsl_vector* K_vi, const NonHumanHostsType& nonHumanHosts)
{
	// \f$P_{dif}\f$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle and gets infected.
  gsl_vector* P_dif;
	// \f$P_{duf}\f$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle and does not get infected.
  gsl_vector* P_duf; 
	
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
  P_dif = gsl_vector_calloc(theta_p);
  P_duf = gsl_vector_calloc(theta_p);


	// We note again, that this code is written assuming there is only
	// one type of hosts. 
	
	// Refer to papers noted above for equations.
        // P_A and P_Ai are described in CalcInitMosqEmergeRate.
  double totalAvailability = alpha_i * N_i;
  /* FIXME
  for (NonHumanHostsType::const_iterator nnh = nonHumanHosts.begin(); nnh != nonHumanHosts.end(); ++nnh)
    totalAvailability += nnh->availability; */
  double P_A = exp(-(totalAvailability + mu_vA)*theta_d);
  double P_Ai = (1-P_A)*(alpha_i*N_i)/(totalAvailability + mu_vA);
	// \f$P_{df}\f$: Probability that a mosquito finds a host on a given
	// night and then completes the feeding cycle.
	// Only part of calculation here; rest after P_df has been used to calculate P_dif and P_duf
  double P_df = P_Ai*P_B_i*P_C_i*P_D_i*P_E_i;

	// Evaluate P_dif and P_duf.
	// We would add P_Ai*P_B_i*P_C_i*P_D_i*P_E_i*kappa for non-human hosts, except kappa is always zero
	
	// P_dif:
  gsl_vector_memcpy(P_dif, K_vi);
  gsl_vector_scale(P_dif, P_df);

	// P_duf:
  gsl_vector_set_all(P_duf, 1.0);
  gsl_vector_sub(P_duf, K_vi);
  gsl_vector_scale(P_duf, P_df);
  /* FIXME
  for (NonHumanHostsType::const_iterator nnh = nonHumanHosts.begin(); nnh != nonHumanHosts.end(); ++nnh)
    P_df += nnh->availability * nnh->probSurviveCycle;
	*/

  sumkplus = 0;
  sumkplusPtr = &sumkplus;

	// Calculate probabilities of mosquito surviving the extrinsic
	// incubation period.]
	// These currently do not depend on the phase of the period.
  CalcPSTS(sumkplusPtr, sumklplus, P_A, P_df);
  sumkplus = *sumkplusPtr;

	/* // Print sumkplus and sumklplus
  trace << "In CalcUpsilonOneHost, sumkplus = " << sumkplus << endl;
  for(l=1; l<=tau-1; l++){
  trace << "sumklplus(" << l << ") = " <<  sumklplus[l-1] << endl;
}
  getchar();
        */ 
        

  size_t indexSv = 2*mt;
	// We start creating the matrices now.
	// Refer to Section 2.1 of JBD Paper for how this matrix is created.
  for (size_t k=0; k < theta_p; k++){
    for (size_t i=0; i < eta; i++){
			// Set 1's along the subdiagonal of all rows except the three
			// rows for the the main system variables.
      if (i!=0uL && i!=mt && i!=indexSv) {
        gsl_matrix_set(Upsilon[k],i,i-1,1.0);
      }
    }

		
		// for \f$N_v\f$.
    gsl_matrix_set(Upsilon[k],0,0,P_A);
    double temp = P_df + gsl_matrix_get(Upsilon[k], 0, tau-1);
    gsl_matrix_set(Upsilon[k],0,tau-1,temp);

		// for \f$O_v\f$.
		// We add theta_p to i, to ensure that it's positive.
		// % is the mod function.
    temp = gsl_vector_get(P_dif,(k+theta_p-tau)%theta_p);		
    gsl_matrix_set(Upsilon[k],mt,tau-1,temp);
    gsl_matrix_set(Upsilon[k],mt,mt,P_A);
    temp = gsl_vector_get(P_duf,(k+theta_p-tau)%theta_p) 
        + gsl_matrix_get(Upsilon[k], mt, mt+tau-1);	
    gsl_matrix_set(Upsilon[k],mt,mt+tau-1,temp);

		// for \f$S_v\f$.
    temp = gsl_vector_get(P_dif,(k+theta_p-theta_s)%theta_p)*sumkplus;
    gsl_matrix_set(Upsilon[k],indexSv,theta_s-1,temp);
    gsl_matrix_set(Upsilon[k],indexSv,mt+theta_s-1,-temp);
    for (size_t l=1; l <= tau-1; l++){
      temp = gsl_vector_get(P_dif,(k+theta_p-theta_s-l)%theta_p)*sumklplus[l-1];
      gsl_matrix_set(Upsilon[k],indexSv, theta_s+l-1, temp);
      gsl_matrix_set(Upsilon[k],indexSv, mt+theta_s+l-1, -temp);
    }
    gsl_matrix_set(Upsilon[k], indexSv, indexSv, P_A);
    temp = P_df + gsl_matrix_get(Upsilon[k], indexSv, indexSv+tau-1);
    gsl_matrix_set(Upsilon[k], indexSv, indexSv+tau-1, temp);


  }

# ifdef VectorTransmission_PRINT_CalcUpsilonOneHost
  // We should try to print some of these matrices out to see what they  look like.
  PrintUpsilon(Upsilon, theta_p, eta, P_A, P_Ai, P_df, P_dif, P_duf);
# endif

	// Reference pointers.
  *PAPtr = P_A;
  *PAiPtr = P_Ai;
	
	// Deallocate memory for vectors
  gsl_vector_free(P_dif);
  gsl_vector_free(P_duf);
}


int CalcSvDiff_rf(const gsl_vector* x, void* p, gsl_vector* f){
	// Cast the incoming pointer, p, to point at a structure of type
	// SvDiffParams.
  SvDiffParams* params = (SvDiffParams *) p;
  VectorEmergence* emerge = params->emerge;
  
  /* If x is the same as lastNv0 (last array used to calculate f), then don't
   * recalculate f. This was intended as an optimisation, but x is at least
   * slightly different each time, so it has no effect. */
  for (size_t i = 0; i < emerge->theta_p; ++i) {
    if (gsl_vector_get(params->lastNv0, i) == gsl_vector_get(x, i))
      continue;	// if same, test next values (or leave for loop)
    // if different, calculate (then break)
    
    // To set f, we simply call CalcSvDiff. It's probably easier than rewriting
    // this code.
    emerge->CalcSvDiff(params->lastS_vDiff, params->S_vFromEIR, x, params->inv1Xtp);
    
    // Calculate the l^1 norm of SvDiff.
    emerge->trace << "The \f$l^1\f$ norm of S_vDiff is " <<  gsl_blas_dasum(params->lastS_vDiff) << endl;
    gsl_vector_memcpy (params->lastNv0, x);
    break;	// we only want to run this code once
  }
  
  emerge->counterSvDiff++;
  emerge->trace << "In CalcSvDiff_rf for the "<< emerge->counterSvDiff <<"th time"<<endl;
  
  gsl_vector_memcpy (f, params->lastS_vDiff);
  
  return GSL_SUCCESS;
}


void VectorEmergence::CalcSvDiff(gsl_vector* S_vDiff, const gsl_vector* S_vFromEIR, 
		const gsl_vector* N_v0, const gsl_matrix* inv1Xtp)
{
	// Periodic orbit of the number of infectious mosquitoes calculated for
	// the given N_v0.
	// \f$S_v\f$.
  gsl_vector* SvfromNv0 = gsl_vector_calloc(theta_p);

	// Calculate the forcing term for each time in the period.
  CalcLambda(N_v0);

	// Calculate the periodic orbit for the given N_v0.
  CalcXP(inv1Xtp);
  
	// Extract the number of infectious mosquitoes from the full periodic
	// orbit.
  size_t indexSv = 2*mt;
  for (size_t i=0; i<theta_p; i++){
    gsl_vector_set(SvfromNv0, i,
                   gsl_vector_get(x_p[i], indexSv));
  }

# ifdef VectorTransmission_PRINT_CalcSvDiff
  PrintVector("SvfromNv0", SvfromNv0);
# endif

	// Subtract S_vFromEIR from SvfromNv0
  gsl_vector_memcpy(S_vDiff,SvfromNv0);
  gsl_vector_sub(S_vDiff, S_vFromEIR);


  gsl_vector_free(SvfromNv0);
}


void VectorEmergence::CalcLambda(const gsl_vector* N_v0)
{
  for(size_t t=0; t < theta_p; t++){
    gsl_vector_set(Lambda[t], 0, gsl_vector_get(N_v0, t));
  }
	
# ifdef VectorTransmission_PRINT_CalcLambda
  // We should try to print some of these vectors out to see what they  look like.
  PrintVector("LambdaFirst", Lambda[0]);
  PrintVector("LambdaLast", Lambda[theta_p-1]);
# endif
}


void VectorEmergence::CalcXP(const gsl_matrix* inv1Xtp)
{
  gsl_vector* vtemp = gsl_vector_calloc(eta);
	// gsl_vector* vtempsum = gsl_vector_calloc(eta);
  gsl_matrix* mtemp = gsl_matrix_calloc(eta, eta);

	// Initial condition for periodic orbit.
  gsl_vector* x0p = gsl_vector_calloc(eta);

	// Evaluate the initial condition of the periodic orbit.
	// Please refer to paper [add reference to paper and equation
	// number here] for the expression for \f$x_0\f$.
  for(size_t i=0; i < theta_p; i++){
    FuncX(mtemp, theta_p, i+1);
    // vtemp = 1.0 * mtemp * Lambda[i] + 1.0 * vtemp :
    gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, const_cast<const gsl_vector*>(Lambda[i]), 1.0, vtemp);
		// gsl_vector_add(vtempsum, vtemp);
  }
  gsl_blas_dgemv(CblasNoTrans, 1.0, inv1Xtp, vtemp, 0.0, x0p);

# ifdef VectorTransmission_PRINT_CalcXP
  PrintVector("x0p", x0p);
# endif

	// We evalute the full periodic orbit now. 
	// Note: to try to keep the indices consistent with our notes and MATLAB, 
	// x_p[0] will refer to x_p(1): because Upsilon[0] refers to Upsilon(1).
	// Thus, x_p[theta_p-1] = x0p. We can check this to make sure.
  for(size_t t=0; t<theta_p; t++){
#ifdef DEBUG_PRINTING
    trace << "t=" << (t+1) << '\r';
#endif
		// gsl_vector_set_zero(vtemp);
		// gsl_vector_set_zero(vtempsum);
    gsl_vector_set_zero (x_p[t]);
    FuncX(mtemp, t+1, 0);
    gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, x0p, 1.0, x_p[t]);
    for(size_t i=0; i<=t; i++){
			// trace << "t=" << t << " i=" <<  i << endl;
      FuncX(mtemp, t+1, i+1);
      gsl_blas_dgemv(CblasNoTrans, 1.0, mtemp, const_cast<const gsl_vector*>(Lambda[i]), 1.0, x_p[t]);
    }
  }
  
# ifdef VectorTransmission_PRINT_CalcXP
  // We should try to print some of these vectors out to see what they  look like.
  PrintXP(x_p, eta, theta_p);
# endif

  gsl_vector_free(vtemp);
	// gsl_vector_free(vtempsum);
  gsl_vector_free(x0p);
  gsl_matrix_free(mtemp);
}


void VectorEmergence::CalcPSTS(double* sumkplusPtr, double* sumklplus, double P_A, double P_df) const
{
  int klplus;	// \f$k_{l+}\f$ in model.
	// klplus = (int *)malloc((tau-1)*sizeof(int)); Define temporarily.
        // \f$k_+\f$ in model:
  int kplus = (int) ((double(theta_s)/tau)-1.); // = floor(theta_s/tau)-1;

	// Evaluate sumkplus
  double sumkplus = 0.;
  for (int j=0; j <= kplus; j++){
    double tempbin = binomial(theta_s-(j+1)*tau+j,j);
    double temppap = pow(P_A,(double)theta_s-(j+1)*tau);
    double temppdfp = pow(P_df,(double)j);
    double temp = tempbin*temppap*temppdfp;
    sumkplus += temp;
  }
  *sumkplusPtr = sumkplus;

	/* // Check results for kplus. 
  itemp = (int) (7./2.) - 1.;
  trace << "kplus = " << kplus << endl;
  trace << "itemp = " << itemp << endl;
  getchar();
        */

	/* // Print sumkplus
  trace << "In CalcPSTS(), sumkplus = " << sumkplus << endl;
  getchar();
        */ 

	// Evaluate sumklplus
  for (size_t l=1; l <= tau-1; l++){
    klplus = (int) (double(theta_s+l)/tau) -2; // = floor((theta_s+l)/tau)-2;
    sumklplus[l-1] = 0;
		// trace << "For l = " << l << ", klplus = " <<  klplus << endl;

    for(int j=0; j<=klplus; j++){
      double tempbin = binomial(theta_s+l-(j+2)*tau+j,j);
      double temppap = pow(P_A,(double)(theta_s+l-(j+2)*tau));
      double temppdfp = pow(P_df,(double)(j+1));
      double temp = tempbin*temppap*temppdfp;
      sumklplus[l-1] += temp;
			// trace << "For j = " << j << ", tempsum = " <<  temp << endl;
    }
		// trace << "sumklplus(" << l << ") = " <<  sumklplus[l-1] << endl;
  }
}


void VectorEmergence::FuncX(gsl_matrix* X, size_t t, size_t s) const
{
  gsl_matrix_set_identity(X);
  
  for (size_t i=s; i<t; i++){
    //Pass the blas function a copy of Upsilon[i] so it can't change the matrix there
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, const_cast<const gsl_matrix*>(Upsilon[i]), X, 0.0, memMatrixEtaSq);
    gsl_matrix_memcpy(X, memMatrixEtaSq);
  }
}


double VectorEmergence::CalcSpectralRadius(const gsl_matrix* A) const
{
  // Copy A into memMatrixEtaSq to keep it safe; gsl_eigen_nonsymm changes its copy
  gsl_matrix_memcpy(memMatrixEtaSq, A);

  // Calculate eigenvalues of memMatrixEtaSq:
  gsl_eigen_nonsymm(memMatrixEtaSq, memVectorComplexEta, memEigenWorkspace);

# ifdef VectorTransmission_PRINT_CalcSpectralRadius
  //logFile << "Eigenvalues =" << endl;
  //gsl_vector_complex_fprintf(fpp, eval, "%e");
# endif

  // Calculate the absolute values of the eigenvalues, and find the largest.
  double sr = 0.0;	// sprectral radius
  for(size_t i=0; i<eta; i++){
    double temp = gsl_complex_abs(gsl_vector_complex_get(memVectorComplexEta, i));
    if (temp > sr)
      sr = temp;
  }
  
  return sr;
}


void VectorEmergence::CalcInv1minusA(gsl_matrix* inv1A, const gsl_matrix* A) const
{
  // We calculate (I-A) in memMatrixEtaSq.
  gsl_matrix_set_identity(memMatrixEtaSq);	// memMatrixEtaSq = I.
  gsl_matrix_sub(memMatrixEtaSq, A);		// memMatrixEtaSq = I-A.

  // Calculate LU decomposition of (I-A), and use to calculate inverse.
  int signum;
  gsl_linalg_LU_decomp(memMatrixEtaSq, memPermutation, &signum);
  gsl_linalg_LU_invert(memMatrixEtaSq, memPermutation, inv1A);

# ifdef VectorTransmission_PRINT_CalcInv1minusA
  char invname[15] = "inv1minusA"; // Name of matrix (when printing to file).
  PrintMatrix(invname, inv1A, eta, eta);
# endif
}


void VectorEmergence::CalSvfromEIRdata(gsl_vector* Sv, double P_Ai, const gsl_vector* Xi_i) const
{
  // Sv(t) = Xi_i(t)*(N_i/(P_Ai*P_B_i))
  gsl_vector_memcpy(Sv, Xi_i);
  gsl_vector_scale(Sv, N_i/(P_Ai*P_B_i));	
}


double VectorEmergence::binomial(int n, int k) const
{
  return gsl_sf_fact((unsigned int) n) / (gsl_sf_fact((unsigned int) k) * gsl_sf_fact((unsigned int) (n-k)));
}



/******************************************************************************
  Printing routines below. Most are only optionally compiled in.
******************************************************************************/

void VectorEmergence::PrintRootFindingStateTS(size_t iter, const gsl_multiroot_fsolver* srootfind) const
{
  trace << "iter = " << (int)iter
	<< "\tN_v0(1) = " << gsl_vector_get(srootfind->x, 0)
	<< "\t||f||_1 = " << gsl_blas_dasum(srootfind->f)	// Calculate the \f$l^1\f$ norm of f.
	<< endl;
}

#ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate	// only use
void VectorEmergence::PrintParameters(size_t theta_p, size_t tau, size_t theta_s,
		double N_i, double alpha_i, double mu_vA,
		double theta_d, double P_B_i, double P_C_i, double P_D_i, double P_E_i,
		const gsl_vector* K_vi, const gsl_vector* Xi_i) const
{
  logFile << "theta_p = " << theta_p << ";" << endl;
  logFile << "tau = " << tau << ";" << endl;
  logFile << "theta_s = " << theta_s << ";" << endl;

  logFile << "N_i = " << N_i << ";" << endl;
  logFile << "alpha_i = " << alpha_i << ";" << endl;
  logFile << "mu_vA = " << mu_vA << ";" << endl;
  logFile << "theta_d = " << theta_d << ";" << endl;
  logFile << "P_B_i = " << P_B_i << ";" << endl;
  logFile << "P_C_i = " << P_C_i << ";" << endl;
  logFile << "P_D_i = " << P_D_i << ";" << endl;
  logFile << "P_E_i = " << P_E_i << ";" << endl;

	
  PrintVector("K_vi", K_vi);

  PrintVector("Xi_i", Xi_i);

  // Let's do this properly.
  for (size_t i=0; i<theta_p; i++){
    logFile << "K_vi(" <<  i+1 << ") = " <<  gsl_vector_get(K_vi, i) << ";" << endl;
  }

  for (size_t i=0; i<theta_p; i++){
    logFile << "Xi_i(" <<  i+1 << ") = " <<  gsl_vector_get(Xi_i,i) << ";" << endl;
  }
}
#endif

#ifdef VectorTransmission_PRINT_CalcUpsilonOneHost	// only use
void VectorEmergence::PrintUpsilon(const gsl_matrix *const * Upsilon, size_t theta_p,
		size_t eta, double P_A, double P_Ai, double P_df,
		const gsl_vector* P_dif, const gsl_vector* P_duf) const
{
  logFile << "P_A = " << P_A << endl;
  logFile << "P_Ai = " << P_Ai << endl;
  logFile << "P_df = " << P_df << endl;

  /*
  for (size_t i=0; i<theta_p; i++)
    logFile << "P_dif(" <<  i+1 << ") = " <<  gsl_vector_get(P_dif, i) << endl;

  for (size_t i=0; i<theta_p; i++)
    logFile << "P_duf(" <<  i+1 << ") = " <<  gsl_vector_get(P_duf, i) << endl;
  */

  // Print some Upsilon[k].
  size_t k=0;

  logFile << "Upsilon[" << k << "] =" << endl;
  for (size_t i=0; i < eta; i++){
    for (size_t j=0; j < eta; j++){
      logFile << gsl_matrix_get(Upsilon[k], i, j) << endl;
    }
    logFile << endl;
  }
  
  k = theta_p - 1;

  logFile << "Upsilon[" << k << "] =" << endl;
  for (size_t i=0; i < eta; i++){
    for (size_t j=0; j < eta; j++){
      logFile << gsl_matrix_get(Upsilon[k], i, j) << endl;
    }
    logFile << endl;
  }
}
#endif

#ifdef VectorTransmission_PRINT_CalcXP	// only use
#include <string.h>

void VectorEmergence::PrintXP(const gsl_vector *const * x_p, size_t eta, size_t theta_p) const
{
  // Print all x_p[t]:
  for(size_t t=0; t<theta_p; t++){
    ostringstream temp;
    temp << "x_p(" << t+1 << ')';
    PrintVector(temp.str().c_str(), x_p[t]);
  }
}
#endif

#if defined VectorTransmission_PRINT_CalcInitMosqEmergeRate || defined VectorTransmission_PRINT_CalcInv1minusA
void VectorEmergence::PrintMatrix(const char matrixname[], const gsl_matrix* A, 
                 size_t RowLength, size_t ColLength) const
{
  logFile << "" << matrixname << " =" << endl;
  for (size_t i=0; i < ColLength; i++){
    for (size_t j=0; j < RowLength; j++){
      logFile << "" << gsl_matrix_get(A, i, j) << endl;
    }
    logFile << endl;
  }
}
#endif

#if defined VectorTransmission_PRINT_CalcInitMosqEmergeRate || defined VectorTransmission_PRINT_CalcSvDiff || defined VectorTransmission_PRINT_CalcXP
void VectorEmergence::PrintVector(const char* vectorname, const gsl_vector* v) const
{
  for (size_t i=0; i < v->size; i++){
    double temp = gsl_vector_get(v, i);
    logFile << vectorname << "(" <<  i+1 << ") = " <<  temp << ";" << endl;
  }
}
#endif

void VectorEmergence::PrintArray(const char* vectorname, const double* v, int n) const {
  for (int i=0; i < n; i++){
    logFile << vectorname << "(" <<  i+1 << ") = " <<  v[i] << ";" << endl;
  }
}
void VectorEmergence::PrintVector(const char* vectorname, const vector<double>& v) const {
  for (unsigned int i=0; i < v.size(); i++){
    logFile << vectorname << "(" <<  i+1 << ") = " <<  v[i] << ";" << endl;
  }
}
