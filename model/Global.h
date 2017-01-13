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

/** The Global module is only for items wanted nearly everywhere.
 * It is not for miscelaneous stuff and should only include headers wanted nearly everywhere.
 */
#ifndef Hmod_Global
#define Hmod_Global

#ifdef _MSC_VER
#pragma warning(disable: 4290)  // disable some warnings on MSVC
// Avoid min/max as macros on windows (breaks numeric_limits<T>::max())
#define NOMINMAX
#endif

// unit32_t and similar, mostly for compatibility with MSVC:
#include <boost/cstdint.hpp>

#include <boost/math/special_functions/fpclassify.hpp>
// this provides: isfinite, isinf, isnan, isnormal
// use these functions like this: (boost::math::isnan)(x)
// extra brackets are necessary since isnan etc. may be defined as macros!
// explicit namespace usage avoids confusion with std versions available sometimes

// Checkpointing and time-step operations are used _everywhere_:
#include "sim.h"    // includes util/checkpoint.h and util/mod.h

// foreach "keyword"
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace OM {

using boost::uint8_t;
using boost::uint16_t;
using boost::uint32_t;
using boost::uint64_t;
using boost::int8_t;
using boost::int16_t;
using boost::int32_t;
using boost::int64_t;

using namespace util::checkpoint;
using util::mod;
using util::mod_nn;

}
#endif
