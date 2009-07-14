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
#include "GSLWrapper.h"
#include "inputData.h"
#include "population.h"

gsl_rng * generator;
const gsl_multimin_fminimizer_type *T;
gsl_multimin_fminimizer *s;
gsl_multimin_function minex_func;


void GSL_SETUP(){
  //use the mersenne twister generator
  generator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (generator, getISeed());
  T = gsl_multimin_fminimizer_nmsimplex;
}

void  GSL_TEARDOWN(){
  gsl_rng_free (generator);
}

double sampleFromLogNormal(double normp, double meanlog, double stdlog){
  // Used for performance reasons. Calling GLS LOG_NORMAL 5 times is 50% slower.

  double zval=W_UGAUSS_PINV(normp);
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

double W_BETA(double a, double b){
  return gsl_ran_beta (generator,a,b);
}

double W_GAUSS(double mean, double std){
  return gsl_ran_gaussian(generator,std)+mean;
}

double  wCalcRSS(const gsl_vector *v, void* params){
  double x, y;  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  return Population::setDemoParameters(x,y); 
}

double w_minimize_calc_rss(double *param1,double *param2){
   
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  const gsl_multimin_fminimizer_type *T =gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_function minex_func;
  double par[2] = { 0.0, 0.0 };
  ss = gsl_vector_alloc (2);

  gsl_vector_set_all (ss, 0.1);

  x = gsl_vector_alloc (2);

  gsl_vector_set (x, 0, *param1);
  gsl_vector_set (x, 1, *param2);
  minex_func.f = &wCalcRSS;
  minex_func.params = (void *)&par;
  minex_func.n = 2;

  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  size_t iter = 0;
  int status;
  double size;
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
      
    if (status) 
      break;

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1e-2);

    if (status == GSL_SUCCESS){
  //do nothing
    }     
  }
  while (status == GSL_CONTINUE && iter < 100);
  double rss = s->fval;
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
  return rss;
}

double W_UGAUSS_P(double x){
  return gsl_cdf_ugaussian_P(x);
}

double W_UGAUSS_PINV(double p){
  return gsl_cdf_ugaussian_Pinv(p);
}

double W_LOGNORMAL(double mean, double std){
  return gsl_ran_lognormal (generator, mean, std);
}

int W_POISSON(double lambda){
   if ((std::fabs(lambda) > std::numeric_limits<double>::max())){
    //This would lead to an inifinite loop in gsl_ran_poisson
     cerr << "lambda isInf" << endl;
     exit(-1);
    }
  return gsl_ran_poisson (generator, lambda);
}

double W_GAMMA(double a, double b){
  /*double xi, eta_m, V2m, V2m_1;
  int m;
  double delta= *a;
  if (delta > 1) {
    xi=gsl_ran_gamma(generator, *a, *b);
  }
  else{    
    m=0;  
    double e=2.71828182;
    double v0 = e/(e + delta);
    do{
      m++;
      V2m_1 = gsl_rng_uniform (generator);
      V2m = gsl_rng_uniform (generator);
      if (V2m_1 <= v0){
        xi = pow(V2m_1,1.0/delta);
        eta_m = V2m*pow(xi,delta - 1);
      }
      else{ 
        xi = 1.0 - log(V2m_1);
        eta_m = V2m*pow(e,-xi);
      };
// TODO: it looks like the following might give more sensible answers if log transformed 
  } while ( eta_m > pow(xi,delta - 1)*pow(e,-xi) );


}
  return xi;
*/
  return gsl_ran_gamma(generator, a, b);
}

double W_UNIFORM(){
  return gsl_rng_uniform (generator);
}

void load_rng_state(int seedFileNumber){
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

void save_rng_state(int seedFileNumber){
  ostringstream seedN;
  seedN << string("seed") << seedFileNumber;
  FILE * f = fopen(seedN.str().c_str(), "wb");
  gsl_rng_fwrite(f, generator);
  fclose (f);
}

