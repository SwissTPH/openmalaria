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


#ifndef GSL_WRAPPER_H
#define GSL_WRAPPER_H

// Note: should be in namespace util too but can't see much point here.
namespace OM {
///@brief A wrapper about some GSL functionality, and some additional functionality based on GSL routines.
namespace gsl {
  ///@brief Random number generators
  //@{
  /** Generate a random number in the range [0,1). */
  double rngUniform();
  
  /** This function returns a Gaussian random variate, with mean mean and standard deviation std. */
  double rngGauss (double mean, double std);
  
  /** This function returns a random variate from the gamma distribution. */
  double rngGamma (double a, double b);
  
  /** This function returns a random variate from the lognormal distribution. */
  double rngLogNormal (double mean, double std);
  
  /** Used for performance reasons. Calling rngLogNormal 5 times is 50% slower. */
  double sampleFromLogNormal (double normp, double meanlog, double stdlog);
  
  /** This function returns a random variate from the beta distribution. */
  double rngBeta(double a, double b);
  
  /** This function returns a random integer from the Poisson distribution with mean lambda. */
  int rngPoisson(double lambda);
  //@}
  
  /** These functions compute the cumulative distribution function P(x) and it's inverse for the unit Gaussian distribution. */
  double cdfUGaussianP (double x);
  /// ditto
  double cdfUGaussianPInv (double p);
  
  /** Presumably from AR. Searches iteratively for a minimum to function func.
   *
   * @param func A function pointer to the function used for ...
   * @param param1 Initial guess of first param taken by func
   * @param param2 Initial guess of second
   */
  void minimizeCalc_rss(double (*func) (double,double), double param1, double param2);
  
  
  ///@brief Setup & cleanup
  //@{
  /// Set up the random-number-generator with seed (usually InputData.getISeed()).
  void setUp(unsigned long int seed);
  /// Free memory
  void tearDown();
  
  void rngSaveState(int seedFileNumber);
  void rngLoadState(int seedFileNumber);
  //@}
}
}
#endif