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


#include "TransmissionModel/VectorInternal.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex_math.h>

const char fntestentopar[30] = "output_ento_para.txt";


void CalcUpsilonOneHost(gsl_matrix** Upsilon, double* PAPtr, 
                        double* PAiPtr, size_t theta_p, size_t eta, size_t mt, size_t tau, 
                        size_t theta_s, size_t n, size_t m, double N_i, double alpha_i, 
                        double mu_vA, double theta_d, double P_B_i, double P_C_i, double P_D_i, 
                        double P_E_i, gsl_vector* K_vi)
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
        

  size_t indexSv = 2*mt;
	// We start creating the matrices now.
	// Refer to Section 2.1 of JBD Paper for how this matrix is created.
  for (size_t k=0; k < theta_p; k++){
    Upsilon[k] = gsl_matrix_calloc(eta, eta);

    for (size_t i=0; i < eta; i++){
			// Set 1's along the subdiagonal of all rows except the three
			// rows for the the main system variables.
      if (i!=0uL && i!=mt && i!=indexSv) {
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
    gsl_matrix_set(Upsilon[k],indexSv,theta_s-1,temp);
    gsl_matrix_set(Upsilon[k],indexSv,mt+theta_s-1,-temp);
    for (size_t l=1; l <= tau-1; l++){
      temp = gsl_vector_get(Pdif,(k+theta_p-theta_s-l)%theta_p)*sumklplus[l-1];
      gsl_matrix_set(Upsilon[k],indexSv, theta_s+l-1, temp);
      gsl_matrix_set(Upsilon[k],indexSv, mt+theta_s+l-1, -temp);
    }
    gsl_matrix_set(Upsilon[k], indexSv, indexSv, P_A);
    temp = Pdf + gsl_matrix_get(Upsilon[k], indexSv, indexSv+tau-1);
    gsl_matrix_set(Upsilon[k], indexSv, indexSv+tau-1, temp);


  }

# ifdef VectorTransmission_PRINT_CalcUpsilonOneHost
  // We should try to print some of these matrices out to see what they  look like.
  PrintUpsilon(Upsilon, theta_p, eta, P_A, P_Ai, Pdf, Pdif, Pduf);
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
  size_t eta = (params->eta);
  size_t mt = (params->mt);
  size_t theta_p = (params->thetap);

	// Recreate a new N_v0 so that we're not restricted by const problems.
  gsl_vector* N_v0 = gsl_vector_calloc(theta_p);
  gsl_vector_memcpy(N_v0, x);


  counterSvDiff++;
  printf("In CalcSvDiff_rf for the %d th time \n", counterSvDiff);

	// To set f, we simply call CalcSvDiff. It's probably easier than rewriting
	// this code.
  CalcSvDiff(f, S_vFromEIR, Upsilon, N_v0, inv1Xtp, 
             eta, mt, theta_p);

	// Calculate the l^1 norm of f.
  SvDiff1norm = gsl_blas_dasum(f);
	
	// gsl_vector_set_all(f, 2.3);

  printf("The $l^1$ norm of S_vDiff is %e \n", SvDiff1norm);
	
  gsl_vector_free(N_v0);
  return GSL_SUCCESS;
}


void CalcSvDiff(gsl_vector* S_vDiff, gsl_vector* S_vFromEIR, 
                gsl_matrix** Upsilon, gsl_vector* N_v0, gsl_matrix* inv1Xtp, 
                size_t eta, size_t mt, size_t theta_p)
{
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
  CalcLambda(Lambda, N_v0, eta, theta_p);

	// Calculate the periodic orbit for the given N_v0.
  CalcXP(x_p, Upsilon, Lambda, inv1Xtp, eta, theta_p);

	// Extract the number of infectious mosquitoes from the full periodic
	// orbit.
  size_t indexSv = 2*mt;
  for (size_t i=0; i<theta_p; i++){
    gsl_vector_set(SvfromNv0, i,
                   gsl_vector_get(x_p[i], indexSv));
  }

# ifdef VectorTransmission_PRINT_CalcSvDiff
  char SvfromNv0name[15] = "SvfromNv0";
  PrintVector(SvfromNv0name, SvfromNv0, theta_p);
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


void CalcLambda(gsl_vector** Lambda, gsl_vector* N_v0, size_t eta,
                size_t theta_p)
{
  for(size_t t=0; t < theta_p; t++){
    Lambda[t] = gsl_vector_calloc(eta);
    gsl_vector_set(Lambda[t], 0, gsl_vector_get(N_v0, t));
  }
	
# ifdef VectorTransmission_PRINT_CalcLambda
  // We should try to print some of these vectors out to see what they  look like.
  PrintLambda(Lambda, eta);		
# endif
}


void CalcXP(gsl_vector** x_p, gsl_matrix** Upsilon, 
            gsl_vector** Lambda, gsl_matrix* inv1Xtp, size_t eta,
            size_t theta_p)
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
# ifdef VectorTransmission_PRINT_CalcXP
  char x0pname[15] = "x0p";
  PrintVector(x0pname, x0p, eta);
# endif

	// We evalute the full periodic orbit now. 
	// Note: to try to keep the indices consistent with our notes and MATLAB, 
	// x_p[0] will refer to x_p(1): because Upsilon[0] refers to Upsilon(1).
	// Thus, x_p[theta_p-1] = x0p. We can check this to make sure.
  for(size_t t=0; t<theta_p; t++){
		// Print t 
    printf("t=%u \r", (unsigned)(t+1));
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

# ifdef VectorTransmission_PRINT_CalcXP
  // We should try to print some of these vectors out to see what they  look like.
  PrintXP(x_p, eta, theta_p);		
# endif

  gsl_vector_free(vtemp);
	// gsl_vector_free(vtempsum);
  gsl_vector_free(x0p);
  gsl_matrix_free(mtemp);
}


void CalcPSTS(double* sumkplusPtr, double* sumklplus, size_t theta_s,
              size_t tau, double P_A, double Pdf)
{
  double taud = (double)tau;
  double thetasd = (double) theta_s;

  int klplus;	// $k_{l+}$ in model.
	// klplus = (int *)malloc((tau-1)*sizeof(int)); Define temporarily.
        // $k_+$ in model:
  int kplus = (int) ((thetasd/taud)-1.); // = floor(theta_s/tau)-1;

	// Evaluate sumkplus
  double sumkplus = 0.;
  for (int j=0; j <= kplus; j++){
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
    klplus = (int) ((thetasd+l)/taud) -2; // = floor((theta_s+l)/tau)-2;
    sumklplus[l-1] = 0;
		// printf("For l = %d, klplus = %d \n", l, klplus);

    for(int j=0; j<=klplus; j++){
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


void FuncX(gsl_matrix* X, gsl_matrix** Upsilon, size_t t, size_t s, size_t eta)
{
  gsl_matrix* temp = gsl_matrix_calloc(eta, eta); 

  gsl_matrix_set_identity(X);

  for (size_t i=s; i<t; i++){
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Upsilon[i], X, 0.0, temp);
    gsl_matrix_memcpy(X, temp);
  }
  gsl_matrix_free(temp);
}


double CalcSpectralRadius(gsl_matrix* A, size_t n)
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

# ifdef VectorTransmission_PRINT_CalcSpectralRadius
  PrintEigenvalues(eval, n);
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


void CalcInv1minusA(gsl_matrix* inv1A, gsl_matrix* A, size_t n)
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


# ifdef VectorTransmission_PRINT_CalcInv1minusA
  char invname[15] = "inv1minusA"; // Name of matrix (when printing to file).
  PrintMatrix(invname, inv1A, n, n);
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
                             size_t theta_p, char fnrootfindingstate[])
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
  printf("iter = %5u N_v0(1) = % .3f ||f||_1 = % .3f \n", (int)iter, Nv0_0, svdiffsum);
  fprintf(fpp, "iter = %5u N_v0(1) = % .3f ||f||_1 = % .3f \n", (int)iter, Nv0_0, svdiffsum);
  fclose(fpp);
}

#ifdef VectorTransmission_PRINT_CalcInitMosqEmergeRate	// only use
void PrintParameters(char fntestentopar[], size_t theta_p, size_t tau, size_t theta_s, 
                    size_t n, size_t m, double N_i, double alpha_i, double mu_vA, 
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
  size_t i;
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

#ifdef VectorTransmission_PRINT_CalcUpsilonOneHost	// only use
void PrintUpsilon(gsl_matrix** Upsilon, size_t theta_p,
                  size_t eta, double P_A, double P_Ai, double Pdf, gsl_vector* Pdif,
                  gsl_vector* Pduf)
{
/* PrintUpsilon() prints the intermediate results while calculating 
  * Upsilon.
  * 
  * All parameters are IN parameters.
 */

  size_t i;
  size_t j;
  size_t k;
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

#ifdef VectorTransmission_PRINT_CalcXP	// only use
#include <string.h>

void PrintXP(gsl_vector** x_p, size_t eta, size_t theta_p)
{
/* PrintXP() prints out values of XP, the periodic orbit.
  * 
  * All parameters are IN parameters.
 */
  size_t t;
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
    PrintVector(xpvecname, x_p[t], eta);
  }

  fclose(fpp);
}
#endif

#ifdef VectorTransmission_PRINT_CalcLambda	// only use
void PrintLambda(gsl_vector** Lambda, size_t eta)
{
/* PrintLambda() prints some values of Lambda.
  * 
  * All parameters are IN parameters.
 */
  size_t t;
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

#ifdef VectorTransmission_PRINT_CalcSpectralRadius	// only use
void PrintEigenvalues(gsl_vector_complex* eval, size_t n)
{
/* PrintEigenvalues() prints eigenvalues to the given file.
  * 
  * All parameters are IN parameters.
 */
	// size_t i;
	// size_t j;
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

#if defined VectorTransmission_PRINT_CalcInitMosqEmergeRate || defined VectorTransmission_PRINT_CalcInv1minusA
void PrintMatrix(char matrixname[], gsl_matrix* A, 
                 size_t RowLength, size_t ColLength)
{

/* PrintMatrix() prints the given matrix to the given file.
  * 
  * All parameters are IN parameters.
 */

  size_t i;
  size_t j;
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

#if defined VectorTransmission_PRINT_CalcInitMosqEmergeRate || defined VectorTransmission_PRINT_CalcSvDiff || defined VectorTransmission_PRINT_CalcXP
void PrintVector(char vectorname[], gsl_vector* v, size_t n)
{
/* PrintVector() prints the given (GSL) vector to the given file.
  * 
  * All parameters are IN parameters.
 */

  FILE* fpp = fopen(fntestentopar, "a");

  for (size_t i=0; i < n; i++){
    double temp = gsl_vector_get(v, i);
    fprintf(fpp, "%s(%d) = %f; \n", vectorname, i+1, temp);
		
  }
	
  fclose(fpp);
}
#endif

void PrintArray(char vectorname[], double* v, int n){
  FILE* fpp = fopen(fntestentopar, "a");
  
  for (int i=0; i < n; i++){
    fprintf(fpp, "%s(%d) = %f; \n", vectorname, i+1, v[i]);
  }
  fclose(fpp);
}
void PrintArray(char vectorname[], vector<double>& v){
  FILE* fpp = fopen(fntestentopar, "a");
  
  for (unsigned int i=0; i < v.size(); i++){
    fprintf(fpp, "%s(%u) = %f; \n", vectorname, i+1, v[i]);
  }
  fclose(fpp);
}
