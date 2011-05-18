/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#ifndef Hmod_TimedDist
#define Hmod_TimedDist

#include "Global.h"
#include "Population.h"
#include "schema/interventions.h"
#include <boost/ptr_container/ptr_vector.hpp>
#include <bitset>

namespace OM {
    using ::boost::ptr_vector;

/** Used to describe which interventions are in use. */
namespace Interventions {
    enum Flags {
        CHANGE_HS,
        CHANGE_EIR,
        VACCINE,        // any vaccine
        MDA,
        IPTI,
        ITN,
        IRS,
        VEC_AVAIL,
        IMMUNE_SUPPRESSION,
        COHORT,
        LARVICIDING,
        R_0_CASE,
        IMPORTED_INFECTIONS,
        UNINFECT_VECTORS,
        SIZE
    };
}

/** Interface of a timed intervention. */
class TimedIntervention {
public:
    /// Create, passing time of deployment
    TimedIntervention(TimeStep deploymentTime);
    virtual ~TimedIntervention() {}
    virtual void deploy () =0;
    TimeStep time;
};

/** Management of interventions deployed on a per-timestep basis. */
class InterventionManager {
public:
    InterventionManager (const scnXml::Interventions& intervElt, OM::Population& population);
    
    /// Returns true if intervention is active
    inline bool isActive( Interventions::Flags intervention ) const{
        assert( intervention < Interventions::SIZE );
        return activeInterventions[intervention];
    }
    
    /** Deploy interventions for this timestep. */
    void deploy ();
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
	nextTimed & stream;
    }
    
private:
    bitset<Interventions::SIZE> activeInterventions;
    // List of all timed interventions. Should be sorted (time weakly increasing).
    ptr_vector<TimedIntervention> timed;
    uint32_t nextTimed;
};

}
#endif
