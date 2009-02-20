/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

/* This file should contain the headers of all routines that we write in the C
 * program.
 */ 

/* We also include library headers here. */ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multiroots.h>
#include "transmissionModel.h"

//! Transmission models, Chitnis et al
class VectorControl :
public TransmissionModel { 

 public:

  VectorControl() {};
  ~VectorControl() {};
	
  //! initialise the main simulation 
  void initMainSimulation (int populationSize); 

  //maxDurIntPhaseEIR is the maximum number of days for which an
  //intervention phase EIR is specified
  static const int ifCalcMosqEmergeRate = 0;	// TODO: Move to XML.
  // This should be 1 for the entomological model code to run.
  // But we can set it to 0 to allow the program to run faster by 
  // skipping EntoModel.cpp.

 private:

  int _populationSize;

  double mosqEmergeRate[daysInYear]; // Note that this is of size daysInYear and not intervalsPerYear.
 
  //! This subroutine converts vectors of length intervalsPerYear to daysInYear. 
  /*! 
    we expect to put the following this into VectorControl class 
 
    Note that ShortArray is assumed to be a pointer to a double array 
    of length intervalsPerYear. We do not explicitly check this. 
 
    For now, we assume that we will use intervalsPerYear and 
    daysInYear as they are defined in global.f. We do not make this 
    subroutine a general subroutine that converts from a given 
    length to another given length. 
 
    TODO: the former documentation said ShortArray should be of type real while FullArray 
    should be of type real*8. 
 
    Note: We also assume that: 
    daysInYear = interval*intervalsPerYear. 
	 
    \param FullArray an array of doubles of length daysInYear 
    \param ShortArray a pointer to a array of doubles of length intervalsPerYear 
    \sa daysInYear, interval, intervalsPerYear 
  */ 
  void convertLengthToFullYear (double FullArray[daysInYear], double* ShortArray); 
	
  //! previously tested Fortran C Communcation 
  /*! 
    This subroutine tests the communication between Fortran and C. 
    We pass an array to C, ask C to add 1 to all the elements and 
    print the new array from Fortran. 
 
    We also pass two matrices to C: A & B. 
    A is a (defined) 2 x 3 matrix 
    B is a (defined) 3 x 2 matrix. 
    In C, we then evaluate 
    C = AB (2 x 2) 
    D = A + B^T (2 x 3) 
 
    We see if C returns to Fortran what we would expect. 
    TODO: when can we scrap this? 
  */ 

  void testFortranCCommunication (); 
 private:
  int i;
  //! get mosquito emergence rates 
  /*! 
    TODO: This function will initializes the system but is not yet implemented. 
    This routine passes the basic entomological parameters (that are already been read, the EIR, and the human infectivity to mosquitoes (all for one type of host) and calculate the mosquito emergenc
    For now, we call testFortranCCommunication to see how well we can pass arrays and matrices between the two. 
    \sa testFortranCCommunication() 
  */ 
  void calMosqEmergeRate (); 

};


//! Formerly tested fortran C interactions, also calculates mosquito emergence rate
/*! 

  testFortranCInteractions() passes test arrays and matrices between C and
  Fortran to see if we can get the communication to work correctly. For now
  we aim to pass an array from Fortran to C, add 1 to all elements of the 
  array, and read the new array in Fortran, while saving a copy of the
  original array in C. 

  We also pass two matrices from Fortran to C:
  A is a (defined) 2 x 3 matrix
  B is a (defined) 3 x 2 matrix.
  In C, we then evaluate:
  C = AB (2 x 2)
  D = A + B^T (2 x 3)
  We see if C returns to Fortran what we would expect.
  TODO: presumably much of this is redundant
  \return The mosquito emergence rate which is as follows
  Units: Mosquitoes/Time
  $N_{v0}$ in model. Vector of length $\theta_p$.
  \sa mosqEmergeRateInitEstimate[]
*/
double TestFortranCInteractions(int* AnnualPeriodPtr, int* nsporePtr,
                                int* FTestArray, int* testArraySizePtr, 
                                double* FAMatrix, int* AMatrixColLengthPtr, int* AMatrixRowLengthPtr, 
                                double* FBMatrix, int* BMatrixColLengthPtr, int* BMatrixRowLengthPtr, 
                                double* FCMatrix, int* CMatrixColLengthPtr, int* CMatrixRowLengthPtr, 
                                double* FDMatrix, int* DMatrixColLengthPtr, int* DMatrixRowLengthPtr);

double CalcInitMosqEmergeRate(double* FMosqEmergeRateVector, int daysInYearPtr,
                              int mosqRestDurationPtr, int EIPDurationPtr, int nHostTypesInitPtr,
                              int nMalHostTypesInitPtr, double popSizeInitPtr, 
                              double hostAvailabilityRateInitPtr, double mosqSeekingDeathRatePtr,
                              double mosqSeekingDurationPtr, double mosqProbBitingPtr,
                              double mosqProbFindRestSitePtr, double mosqProbRestingPtr,
                              double mosqProbOvipositingPtr, double* FHumanInfectivityInitVector,
                              double* FEIRInitVector, double* FMosqEmergeRateInitEstimateVector,
                              char fnametestentopar[]);

void CalcUpsilonOneHost(gsl_matrix** Upsilon, double* PAPtr, 
                        double* PAiPtr, int thetap, int eta, int mt, int tau, 
                        int thetas, int n, int m, double Ni, double alphai, 
                        double muvA, double thetad, double PBi, double PCi, double PDi, 
                        double PEi, gsl_vector* Kvi, char fntestentopar[]);

int CalcSvDiff_rf(const gsl_vector* x, void* p, gsl_vector* f);

// NOTE: Unused
//static int staticCalcSvDiff_rf(const gsl_vector* x, void* p, gsl_vector* f,void* obj);

void CalcSvDiff(gsl_vector* SvDiff, gsl_vector* SvfromEIR, 
                gsl_matrix** Upsilon, gsl_vector* Nv0, gsl_matrix* inv1Xtp, 
                int eta, int mt, int thetap, char fntestentopar[]);

void CalcLambda(gsl_vector** Lambda, gsl_vector* Nv0, int eta,
                int thetap, char fntestentopar[]);

void CalcXP(gsl_vector** xp, gsl_matrix** Upsilon, 
            gsl_vector** Lambda, gsl_matrix* inv1Xtp, int eta,
            int thetap, char fntestentopar[]);

void CalcPSTS(double* sumkplusPtr, double* sumklplus, int thetas,
              int tau, double PA, double Pdf);

void FuncX(gsl_matrix* X, gsl_matrix** Upsilon, int t, int s, int n);

double CalcSpectralRadius(gsl_matrix* A, int n, char fntestentopar[]);

void CalcInv1minusA(gsl_matrix* inv1A, gsl_matrix* A, int n, char fntestentopar[]);

void CalSvfromEIRdata(gsl_vector* Sv, double PAi, double PBi, double Ni, 
                      gsl_vector* Xii);

double binomial(int n, int k);


void CalcCGSLVectorFromCArray(gsl_vector* CGSLVector, double* CArray, 
                              int Length);

void CalcCArrayFromCGSLVector(gsl_vector* CGSLVector, double* CArray, 
                              int Length);

void CalcCGSLMatrixFromFortranArray(gsl_matrix* CMatrix, double* FArray, 
                                    int ColLength, int RowLength);

void CalcFortranArrayFromCGSLMatrix(gsl_matrix* CMatrix, double* FArray, 
                                    int ColLength, int RowLength);

void CalcCGSLVectorFromFortranArray(gsl_vector* CVector, double* FArray, 
                                    int Length);

void CalcFortranArrayFromCGSLVector(gsl_vector* CVector, double* FArray, 
                                    int Length);

void PrintRootFindingStateTS(size_t iter, gsl_multiroot_fsolver* srootfind, 
                             int thetap, char fnrootfindingstate[]);

void PrintParameters(char fntestentopar[], int thetap, int tau, int thetas, 
                     int n, int m, double Ni, double alphai, double muvA, 
                     double thetad, double PBi, double PCi, double PDi, double PEi, 
                     gsl_vector* Kvi, gsl_vector* Xii, gsl_vector* Nv0guess);

void PrintUpsilon(char fntestentopar[], gsl_matrix** Upsilon, int thetap,
                  int eta, double PA, double PAi, double Pdf, gsl_vector* Pdif,
                  gsl_vector* Pduf);

void PrintXP(gsl_vector** xp, int eta, int thetap, char fntestentopar[]);

void PrintLambda(gsl_vector** Lambda, int eta, char fntestentopar[]);

void PrintEigenvalues(char fntestentopar[], gsl_vector_complex* eval, int n);

void PrintMatrix(char fntestentopar[], char matrixname[], gsl_matrix* A, 
                 int RowLength, int ColLength);


void PrintVector(char fntestentopar[], char vectorname[], gsl_vector* v, int n);
