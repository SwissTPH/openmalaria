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

#ifndef Hmod_CommonInfection
#define Hmod_CommonInfection

#include "WithinHost/Infection/Infection.h"

namespace OM { namespace WithinHost {

/** Represent infections used by CommonWithinHost.
 * 
 * All these use a 1-day time step, however CommonWithinHost handles the
 * conversion when the main simulation uses a 5 day time step by updating
 * infections and the PK-PD model multiple times.
 * 
 * Note therefore that sim::ts0(), ts1(), etc. may not always be accurate since
 * it is only updated once per main time step. */
class CommonInfection : public Infection {
public:
    /// @brief Construction and destruction
    //@{
    /// For checkpointing (don't use for anything else)
    CommonInfection(istream& stream);
    /// Per instance initialisation; create new inf.
    CommonInfection(uint32_t genotype) :
	m_genotype(genotype)
    {}
    virtual ~CommonInfection() {}
    //@}
    
    /** Get the infection's genotype. */
    uint32_t genotype()const{ return m_genotype; }
    
    /** Update: calculate new density. Call this once per day.
     * 
     * @param survivalFactor Density multiplier to introduce drug & vaccine effects
     * @param now The simulation time. Use this instead of sim::ts1().
     * @param body_mass Body mass of host in kg
     * @returns True when the infection goes extinct. */
    inline bool update( double survivalFactor, SimTime now, double body_mass ){
	SimTime bsAge = now - m_startDate - latentP;	// age of post-latent-period blood stage
	if( bsAge < sim::zero() )
	    return false;	// latent period (liver stage) — don't do anything
	else
	    return updateDensity( survivalFactor, bsAge, body_mass );
    }
    
protected:
    /** Update: calculate new density.
     *
     * @param survivalFactor Density multiplier to introduce drug & vaccine
     *   effects
     * @param bsAge Age of the patent blood-stage infection (sim::zero() on
     *  first day). Note that liver and pre-patent blood stages occur before
     *  this, but this function is not called during those stages.
     * @param body_mass Body mass of host in kg
     * @returns True when the infection goes extinct. */
    virtual bool updateDensity( double survivalFactor, SimTime bsAge, double body_mass )=0;
    
    virtual void checkpoint (ostream& stream);
    
private:
    /// Genotype of infection (a code; see Genotypes class).
    uint32_t m_genotype;
};

} }
#endif