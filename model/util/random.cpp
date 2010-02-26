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

/* This module contains the random-number generator and distributions wrapper.
 *
 * Currently the implementation uses GSL, though I did look into migrating it to
 * boost. The current boost code is only an addition to the GSL code, however.
 */

#define OM_RANDOM_USE_BOOST

#include "util/random.h"
#include "util/errors.hpp"
#include "Global.h"

#ifdef OM_RANDOM_USE_BOOST
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/static_assert.hpp>
#endif

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <cmath>
#include <sstream>

// Note: since we're using both gsl and boost files, we should be careful to
// avoid name conflicts. So probably don't use "using namespace boost;".


namespace OM { namespace util {

# ifdef OM_RANDOM_USE_BOOST
    static boost::mt19937 boost_generator;
    static boost::uniform_01<boost::mt19937&> rng_uniform01 (boost_generator);
    
    long unsigned int boost_rng_get (void*) {
	BOOST_STATIC_ASSERT (sizeof(uint32_t) <= sizeof(long unsigned int));
	return static_cast<long unsigned int> (boost_generator ());
    }
    double boost_rng_get_double_01 (void*) {
	return rng_uniform01 ();
    }
    
    static const gsl_rng_type boost_mt_type = {
	"boost_mt19937",		// name
	boost_generator.max(),	// max value
	boost_generator.min(),	// min value
	0,					// size of state; not used here
	NULL,				// re-seed function; don't use
	&boost_rng_get,
	&boost_rng_get_double_01
    };
# endif

// This should be created and deleted automatically, taking care of
// allocating and freeing the generator.
struct generator_factory {
    gsl_rng * gsl_generator;
    
    generator_factory () {
#	ifdef OM_RANDOM_USE_BOOST
	// In this case, I construct a wrapper around boost's generator. The reason for this is
	// that it allows use of distributions from both boost and GSL.
	gsl_generator = new gsl_rng;
	gsl_generator->type = &boost_mt_type;
	gsl_generator->state = NULL;
#	else
	//use the mersenne twister generator
	gsl_generator = gsl_rng_alloc(gsl_rng_mt19937);
#	endif
    }
    ~generator_factory () {
#	ifdef OM_RANDOM_USE_BOOST
	delete gsl_generator;
#	else
	gsl_rng_free (gsl_generator);
#	endif
    }
} rng;

// -----  set-up, tear-down and checkpointing  -----

void random::seed (uint32_t seed) {
# ifdef OM_RANDOM_USE_BOOST
    if (seed == 0) seed = 4357;	// gsl compatibility − ugh
    boost_generator.seed (seed);
# else
    gsl_rng_set (rng.gsl_generator, seed);
# endif
}

void random::checkpoint (istream& stream, int seedFileNumber) {
# ifdef OM_RANDOM_USE_BOOST
    // Don't use OM::util::checkpoint function for loading a stream; checkpoint::validateListSize uses too small a number.
    string str;
    size_t len;
    len & stream;
    str.resize (len);
    stream.read (&str[0], str.length());
    if (!stream || stream.gcount() != streamsize(len))
	throw checkpoint_error ("stream read error string");
    istringstream ss (str);
    ss >> boost_generator;
# else
    
    ostringstream seedN;
    seedN << string("seed") << seedFileNumber;
    FILE * f = fopen(seedN.str().c_str(), "rb");
    if (f == NULL)
	throw checkpoint_error (string("load_rng_state: file not found: ").append(seedN.str()));
    if (gsl_rng_fread(f, rng.gsl_generator) != 0)
	throw checkpoint_error ("gsl_rng_fread failed");
    fclose (f);
# endif
}

void random::checkpoint (ostream& stream, int seedFileNumber) {
# ifdef OM_RANDOM_USE_BOOST
    ostringstream ss;
    ss << boost_generator;
    ss.str() & stream;
# else
    
    ostringstream seedN;
    seedN << string("seed") << seedFileNumber;
    FILE * f = fopen(seedN.str().c_str(), "wb");
    if (gsl_rng_fwrite(f, rng.gsl_generator) != 0)
	throw checkpoint_error ("gsl_rng_fwrite failed");
    fclose (f);
# endif
}


// -----  random number generation  -----

double random::uniform_01 () {
    //TODO: check both GSL and boost produce the same result here
    return gsl_rng_uniform (rng.gsl_generator);
    //return rng.rng_uniform01 ();
}

double random::gauss (double mean, double std){
    return gsl_ran_gaussian(rng.gsl_generator,std)+mean;
}

double random::gamma (double a, double b){
    return gsl_ran_gamma(rng.gsl_generator, a, b);
}

double random::log_normal (double mean, double std){
    return gsl_ran_lognormal (rng.gsl_generator, mean, std);
}

double random::sampleFromLogNormal(double normp, double meanlog, double stdlog){
    // Used for performance reasons. Calling GLS LOG_NORMAL 5 times is 50% slower.
    
    double zval = gsl_cdf_ugaussian_Pinv (normp);
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

double random::beta (double a, double b){
    return gsl_ran_beta (rng.gsl_generator,a,b);
}

int random::poisson(double lambda){
    if (!finite(lambda)) {
	//This would lead to an inifinite loop in gsl_ran_poisson
	cerr << "lambda isInf" << endl;
	exit(-1);
    }
    return gsl_ran_poisson (rng.gsl_generator, lambda);
}

} }
