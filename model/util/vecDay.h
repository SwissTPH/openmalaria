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

#ifndef Hmod_OM_util_vecDay
#define Hmod_OM_util_vecDay

#include "Global.h"
#include <vector>

// #if __cplusplus >= 201402L
// // Update classes to mimic std::vector members
// #error "vecDay should be updated if switching to a later C++ version"
// #endif

namespace OM {
namespace util {

namespace vectors{

  /** The inverse of logDFT (or an approximation, when N&lt;T or
   * tArray.size() ≠ T). Result may also be rotated.
   * 
   * (This was called calcExpFourierSeries, and does essentially the same
   * thing.)
   *
   * @param tArray Array to fill with exponated values from Fourier series.
   * Length should already be set. Need not have the same length as the
   * array used to calculate FC.
   * @param FC Fourier coefficients (a0, a1,b1, a2,b2, ...); can be any length
   * so long as it is odd.
   * @param rAngle Angle to rotate generated series by in radians: [0,2π] */
  void expIDFT (std::vector<double>& tArray, const vector<double>& FC, double rAngle);
}

}
}
#endif
