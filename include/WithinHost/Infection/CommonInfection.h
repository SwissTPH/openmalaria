/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_CommonInfection
#define Hmod_CommonInfection

#include "WithinHost/Infection/Infection.h"

namespace OM { namespace WithinHost {

/** Represent infections used by CommonWithinHost.
 * 
 * All these use a 1-day time-step, however CommonWithinHost also allows them
 * to be used where the simulation timestep is set to 5 days, therefore
 * implementors should not assume TimeStep::interval is one and should be
 * careful when using TimeStep::simulation. */
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
    
    /** Update: calculate new density. Call this once per day.
     * 
     * @param survivalFactor Density multiplier to introduce drug & vaccine effects
     * @param day Where TimeStep::interval is 1, this should be 0. Where a
     *  timestep is more than one day long (and thus this function is called
     *  multiple times per timestep), day should be the call number, 0, 1, ...
     *  up to TimeStep::interval - 1.
     * @returns True when the infection goes extinct. */
    inline bool update (double survivalFactor, int day){
	TimeStep ageTS = TimeStep::simulation - _startdate - latentp;	// age of post-latent-period blood stage
	if( ageTS < TimeStep(0) )
	    return false;	// latent period (liver stage) — don't do anything
	else
	    return updateDensity( survivalFactor, ageTS.inDays() + day );
    }
    
protected:
    /** Update: calculate new density.
     *
     * @param survivalFactor Density multiplier to introduce drug & vaccine
     *   effects
     * @param ageDays Age of the patent blood-stage infection in
     *   days (0 on first day). Note that liver and pre-patent blood stages
     *   occur before this, but this function is not called during those stages.
     * @returns True when the infection goes extinct. */
    virtual bool updateDensity (double survivalFactor, int ageDays) =0;
};

} }
#endif