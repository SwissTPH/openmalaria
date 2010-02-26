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
namespace rng {
    ///@brief Setup & cleanup; checkpointing
    //@{
    /// Reseed the random-number-generator with seed (usually InputData.getISeed()).
    void seed (uint32_t seed);
    
    //TODO: checkpoint
    //@}
    
    ///@brief Random number distributions
    //@{
    /** Generate a random number in the range [0,1). */
    double uniform01 ();
    //@}
};
}
// convenience: make rng available in the ::OM namespace
namespace rng = util::rng;
}
