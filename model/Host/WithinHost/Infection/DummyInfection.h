/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
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

#ifndef Hmod_DummyInfection
#define Hmod_DummyInfection

#include "Host/WithinHost/Infection/CommonInfection.h"

namespace OM { namespace WithinHost {

//! Create a dummy infection (used directly by unit test).
CommonInfection* createDummyInfection (LocalRng& rng, uint32_t protID, int origin);

//!  Models of infection.
/*!
  Models related to the within-host dynamics of infections.
*/
class DummyInfection : public CommonInfection {
public:
    /// For checkpointing (don't use for anything else)
    DummyInfection (istream& stream);
    //! Constructor
    DummyInfection (LocalRng& rng, uint32_t protID, int origin);
    
    virtual ~DummyInfection () {}
    
    static void init ();
    
    virtual bool updateDensity( LocalRng& rng, double survivalFactor, SimTime bsAge, double );
};

} }
#endif
