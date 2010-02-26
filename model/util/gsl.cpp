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

#include "util/gsl.h"
#include "inputData.h"
#include "Global.h"
#include "util/errors.hpp"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
using namespace std;

namespace OM { namespace util {


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

} }