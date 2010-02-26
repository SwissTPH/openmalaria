/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "util/random.h"
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <sstream>

using namespace boost;

namespace OM { namespace util {

namespace random {
    mt19937 generator;
    uniform_real<> dist_uniform01 (0,1);
    variate_generator<mt19937&, uniform_real<> > rng_uniform01 (generator, dist_uniform01);
    
    
    // -----  set-up, tear-down and checkpointing  -----
    
    void seed (uint32_t seed) {
	generator.seed (seed);
    }
    
    void checkpoint (istream& stream) {
	string str;
	str & stream;
	istringstream ss (str);
	ss >> generator;
    }
    void checkpoint (ostream& stream) {
	ostringstream ss;
	ss << generator;
	ss.str() & stream;
    }
    
    // -----  random number generation  -----
    
    double uniform01 () {
	return rng_uniform01 ();
    }
    /*
    double Random::rngGauss (double mean, double std){
	return gsl_ran_gaussian(generator,std)+mean;
    }
    
    double Random::rngGamma (double a, double b){
	return gsl_ran_gamma(generator, a, b);
    }
    
    double Random::rngLogNormal (double mean, double std){
	return gsl_ran_lognormal (generator, mean, std);
    }
    
    double Random::sampleFromLogNormal(double normp, double meanlog, double stdlog){
	// Used for performance reasons. Calling GLS LOG_NORMAL 5 times is 50% slower.
	
	double zval = Random::cdfUGaussianPInv (normp);
	/*
	Why not zval=W_UGAUSS?
	where normp is distributed uniformly over [0,1],
	zval is distributed like a standard normal distribution
	where normp has been transformed by raising to the power of 1/(T-1) 
	zval is distributed like a uniform gauss	times 4* F(x,0,1)^3, where F(x,0,1) ist the cummulative
	distr. function of a uniform gauss
	* /
	return exp(meanlog+stdlog*((float)zval));
    }
    
    double Random::rngBeta (double a, double b){
	return gsl_ran_beta (generator,a,b);
    }
    
    int Random::rngPoisson(double lambda){
	if (!finite(lambda)) {
	    //This would lead to an inifinite loop in gsl_ran_poisson
	    cerr << "lambda isInf" << endl;
	    exit(-1);
	}
	return gsl_ran_poisson (generator, lambda);
    }
    
    
    // -----  ?  -----
    
    double Random::cdfUGaussianP (double x){
	return gsl_cdf_ugaussian_P(x);
    }
    
    double Random::cdfUGaussianPInv (double p){
	return gsl_cdf_ugaussian_Pinv(p);
    }
    */
    
    // -----  Setup & cleanup  -----
    /*
    void Random::rngLoadState (int seedFileNumber){
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
    
    void Random::rngSaveState (int seedFileNumber){
	ostringstream seedN;
	seedN << string("seed") << seedFileNumber;
	FILE * f = fopen(seedN.str().c_str(), "wb");
	gsl_rng_fwrite(f, generator);
	fclose (f);
    }*/
};
} }
