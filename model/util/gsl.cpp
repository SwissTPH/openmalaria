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

#include <string>
using namespace std;
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "util/gsl.h"
#include "inputData.h"
#include "Global.h"
#include "util/errors.hpp"

namespace OM { 

gsl_rng * generator;


// -----  random number generation  -----

double gsl::rngUniform (){
  return gsl_rng_uniform (generator);
}

double gsl::rngGauss (double mean, double std){
  return gsl_ran_gaussian(generator,std)+mean;
}

double gsl::rngGamma (double a, double b){
  return gsl_ran_gamma(generator, a, b);
}

double gsl::rngLogNormal (double mean, double std){
  return gsl_ran_lognormal (generator, mean, std);
}

double gsl::sampleFromLogNormal(double normp, double meanlog, double stdlog){
  // Used for performance reasons. Calling GLS LOG_NORMAL 5 times is 50% slower.
  
  double zval = gsl::cdfUGaussianPInv (normp);
    /*
  Why not zval=W_UGAUSS?
  where normp is distributed uniformly over [0,1],
  zval is distributed like a standard normal distribution
  where normp has been transformed by raising to the power of 1/(T-1) 
  zval is distributed like a uniform gauss	times 4* F(x,0,1)^3, where F(x,0,1) ist the cummulative
  distr. function of a uniform gauss
    */
  return exp(meanlog+stdlog*((float)zval));
}

double gsl::rngBeta (double a, double b){
  return gsl_ran_beta (generator,a,b);
}

int gsl::rngPoisson(double lambda){
  if (!finite(lambda)) {
    //This would lead to an inifinite loop in gsl_ran_poisson
    cerr << "lambda isInf" << endl;
    exit(-1);
  }
  return gsl_ran_poisson (generator, lambda);
}


// -----  ?  -----

double gsl::cdfUGaussianP (double x){
  return gsl_cdf_ugaussian_P(x);
}

double gsl::cdfUGaussianPInv (double p){
  return gsl_cdf_ugaussian_Pinv(p);
}

// function pointer for wCalcRSS's func
double (*wCalcRSSFunc) (double param1, double param2);
double wCalcRSS(const gsl_vector *v, void* params){
  double x, y;  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  return (*wCalcRSSFunc)(x,y); 
}
// Not thread safe:
void gsl::minimizeCalc_rss(double (*func) (double,double), double param1,double param2){
    wCalcRSSFunc = func;
    
    gsl_vector *stepSize = gsl_vector_alloc (2);
    gsl_vector_set_all (stepSize, 0.1);
    
    gsl_vector *initialValues = gsl_vector_alloc (2);
    gsl_vector_set (initialValues, 0, param1);
    gsl_vector_set (initialValues, 1, param2);
    
    gsl_multimin_function minex_func;
    minex_func.f = &wCalcRSS;
    minex_func.params = (void *) NULL;
    minex_func.n = 2;
    
    const gsl_multimin_fminimizer_type *T =gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *minimizer = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (minimizer, &minex_func, initialValues, stepSize);
    
    for (size_t iter = 0; iter < 100; ++iter) {
	if (gsl_multimin_fminimizer_iterate(minimizer))
	    throw runtime_error ("gsl_multimin_fminimizer_iterate failed");
	
	double size = gsl_multimin_fminimizer_size (minimizer);
	int status = gsl_multimin_test_size (size, 1e-2);
	if (status == GSL_SUCCESS)
	    break;
    }
    // Call again to set final value. NOTE: this changes previous results.
    wCalcRSS (minimizer->x, NULL);
    
    gsl_vector_free(initialValues);
    gsl_vector_free(stepSize);
    gsl_multimin_fminimizer_free (minimizer);
}


// -----  Setup & cleanup  -----

void gsl::setUp (){
  //use the mersenne twister generator
  generator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (generator, InputData.getISeed());
}

void gsl::tearDown (){
  gsl_rng_free (generator);
}

void gsl::rngLoadState (int seedFileNumber){
  ostringstream seedN;
  seedN << string("seed") << seedFileNumber;
  FILE * f = fopen(seedN.str().c_str(), "rb");
  if (f!=NULL){
    gsl_rng_fread(f, generator);
    fclose (f);
  }else{
    throw runtime_error (string("load_rng_state: file not found: ").append(seedN.str()));
  }
}

void gsl::rngSaveState (int seedFileNumber){
  ostringstream seedN;
  seedN << string("seed") << seedFileNumber;
  FILE * f = fopen(seedN.str().c_str(), "wb");
  gsl_rng_fwrite(f, generator);
  fclose (f);
}

}