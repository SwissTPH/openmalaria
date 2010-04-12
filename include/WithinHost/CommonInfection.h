/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef Hmod_CommonInfection
#define Hmod_CommonInfection

#include "WithinHost/Infection.h"

namespace OM { namespace WithinHost {
    
class CommonInfection : public Infection {
public:
    /// @brief Construction and destruction
    //@{
    /// For checkpointing (don't use for anything else)
    CommonInfection(istream& stream) :
	Infection(stream)
    {}
    /// Per instance initialisation; create new inf.
    CommonInfection(uint32_t protID) :
	Infection(protID)
    {}
    virtual ~CommonInfection() {}
    //@}
    
    //! Get the density of the infection
    inline double getDensity() { return _density; };
    
    /** Update: calculate new density.
    *
    * Currently sets _laggedLogDensities[0] to a large negative number when the
    * infection goes extinct.
    * 
    * @param simulationTime Simulation timestep (expected to be a 1-day timestep)
    * @param survivalFactor Density multiplier to introduce drug & vaccine effects
    * @returns True when the infection goes extinct. */
    virtual bool updateDensity (int simulationTime, double survivalFactor) =0;
};

} }
#endif