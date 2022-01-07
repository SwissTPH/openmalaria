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

#include "util/vectors.h"
#include "util/errors.h"
#include <cstring>
#include <cmath>
#include <cassert>

namespace OM { namespace util {

/*double vectors::sum (const gsl_vector *vec) {
  double r = 0.0;
  for(size_t i = 0; i < vec->size; ++i)
    r += gsl_vector_get( vec, i );
  return r;
}
*/


// NOTE: could replace these algorithms with FFTs, but since these aren't
// called in performance-critical code, the difference would be insignificant.





} }
