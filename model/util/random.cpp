/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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
 * Currently both the GSL and boost generators are implemented. The
 * distributions all come from the GSL library so far.
 * 
 * Using the boost generator appears (in rough tests) to be slightly
 * slower, which is understandable since the GSL distributions must then use a
 * wrapper around the boost generator.
 * 
 * Note: using boost distributions elsewhere could ideally be implemented a
 * little differently, since the distribution objects could in many cases last
 * the length of the program rather than be created on each use.
 */

// Define to use boost distributions. Unfortunately these are not necessarily
// value-stable or portable.
// #define OM_RANDOM_USE_BOOST_DIST

#include "util/random.h"

// Note: since we're using both gsl and boost files, we should be careful to
// avoid name conflicts. So probably don't use "using namespace boost;".


namespace OM { namespace util {

RNG<pcg32> global_RNG;
RNG<ChaCha<8>> master_RNG;

} }
