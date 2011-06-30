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

#include "WithinHost/Infection/Infection.h"

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
    CommonInfection(TimeStep now, uint32_t protID) :
	Infection(now, protID)
    {}
    virtual ~CommonInfection() {}
    //@}
    
    /** Update: calculate new density.
    *
    * @param survivalFactor Density multiplier to introduce drug & vaccine effects
    * @returns True when the infection goes extinct. */
    inline bool update (double survivalFactor){
	TimeStep ageOfInfection = TimeStep::simulation1() - _startdate - latentp;	// age in days
	if( ageOfInfection < TimeStep(0) )
	    return false;	// latent period (liver stage) — don't do anything
	else
	    return updateDensity( survivalFactor, ageOfInfection );
    }
    
protected:
    /** Update: calculate new density.
    *
    * @param survivalFactor Density multiplier to introduce drug & vaccine effects
    * @param ageOfInfection Age of the blood-stage infection in timesteps (0 on first day)
    * @returns True when the infection goes extinct. */
    virtual bool updateDensity (double survivalFactor, TimeStep ageOfInfection) =0;
};

} }
#endif