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

#include "Global.h"

namespace OM { namespace util {

/** Random number generator.
 *
 * This interface should be independant of implementation. */
namespace random {
    ///@brief Setup & cleanup; checkpointing
    //@{
    /// Reseed the random-number-generator with seed (usually InputData.getISeed()).
    void seed (uint32_t seed);
    
    void checkpoint (istream& stream, int seedFileNumber);
    void checkpoint (ostream& stream, int seedFileNumber);
    //@}
    
    ///@brief Random number distributions
    //@{
    /** Generate a random number in the range [0,1). */
    double uniform_01 ();
    
    /** This function returns a Gaussian random variate, with mean mean and standard deviation std. */
    double gauss (double mean, double std);
    
    /** This function returns a random variate from the gamma distribution. */
    double gamma (double a, double b);
    
    /** This function returns a random variate from the lognormal distribution. */
    double log_normal (double mean, double std);
    
    /** Used for performance reasons. Calling rngLogNormal 5 times is 50% slower. */
    double sampleFromLogNormal (double normp, double meanlog, double stdlog);
    
    /** This function returns a random variate from the beta distribution. */
    double beta(double a, double b);
    
    /** This function returns a random integer from the Poisson distribution with mean lambda. */
    int poisson(double lambda);
    //@}
}
} }
