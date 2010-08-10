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
#ifndef Hmod_HostCtsDist
#define Hmod_HostCtsDist

#include "Global.h"

namespace OM { namespace Host {

class Human;	// we have a circular dependency...

/** */
class ContinuousIntervention {
public:
    /** Read XML descriptions.
     *
     * Member-function pointers are used to deploy interventions; these must
     * currently be passed from the Human class. */
    static void init (
	void (Human::*deployItn) (),
	void (Human::*deployIpti) (),
	void (Human::*deployCohort) ()
    );
    
    ContinuousIntervention () :
	nextCtsDist(0)
    {}
    
    /// Deploy any interventions intended for this age in timesteps.
    void deploy (Human* human, int ageTimesteps);
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
	nextCtsDist & stream;
    }
    
private:
    struct AgeIntervention {
	int ageTimesteps;
	double coverage;
	// Member function pointer to the function (in Human) responsible for deploying intervention:
	void (Human::*deploy) ();
	
	inline bool operator< (const AgeIntervention& that) const{
	    return this->ageTimesteps < that.ageTimesteps;
	}
    };
    static vector<AgeIntervention> ctsIntervs;
    
    uint32_t nextCtsDist;
};

} }
#endif
