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
#ifndef Hmod_interventions_HumanComponents
#define Hmod_interventions_HumanComponents

#include "Global.h"
#include "interventions/Interfaces.hpp"
#include "util/DecayFunction.h"

namespace OM {
namespace Host {
    class Human;
}
namespace interventions {
    class VaccineEffect;

// ———  cohorts  ———

namespace Cohort{
    /** Identifiers for types of cohort removal option, followed by NUM in last
     * place (as a counter).
     * 
     * REMOVE_AT_FIRST_BOUT: remove the human from the cohort at the start of
     * each clinical event, if present in the cohort.
     * 
     * REMOVE_AT_FIRST_INFECTION: remove patent humans from the cohort during
     * each survey.
     * 
     * REMOVE_AT_FIRST_TREATMENT: remove the human from the cohort when
     * assigning each course of treatment, if present in the cohort.
     */
    enum RemoveAtCode {
        REMOVE_AT_FIRST_BOUT,
        REMOVE_AT_FIRST_INFECTION,
        REMOVE_AT_FIRST_TREATMENT,
        REMOVE_AT_NUM };
}

// ———  vaccines  ———

namespace Vaccine{
    enum Types { PEV, BSV, TBV, NumVaccineTypes };
}

/** Per vaccine effect (type), per human details. */
class PerEffectPerHumanVaccine {
public:
    //Note: this constructor is only for checkpointing
    PerEffectPerHumanVaccine();
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        effect & stream;
        numDosesAdministered & stream;
        timeLastDeployment & stream;
        initialEfficacy & stream;
        hetSample & stream;
    }
    
private:
    PerEffectPerHumanVaccine( EffectId id, const VaccineEffect& params );
    
    EffectId effect;      // id of effect (for finding parameters)
    // Type of vaccine
    /** Number of vaccine doses this individual has received.
     *
     * If an individual misses one EPI (continuous) vaccine dose, it's
     * intentional that they also miss following EPI doses (unless a timed mass
     * vaccination reintroduces them to the EPI schedule). */
    uint32_t numDosesAdministered;
    /// Timestep of last vaccination with this vaccine type
    TimeStep timeLastDeployment;
    /// Efficacy at last deployment (undecayed)
    double initialEfficacy;
    util::DecayFuncHet hetSample;
    
    friend class PerHumanVaccine;
};

/** Per-human vaccine code. */
class PerHumanVaccine {
public:
    PerHumanVaccine() {}
    
    /** Get one minus the efficacy of the vaccine (1 for no effect, 0 for full effect). */
    double getFactor( Vaccine::Types type )const;
    
    /** Vaccinate unless the passed VaccineLimits specify not to.
     * 
     * @returns true when the vaccine is administered */
    bool possiblyVaccinate( const Host::Human& human, EffectId effectId,
                           interventions::VaccineLimits vaccLimits );

#if 0
    /// Hack for R_0 experiment: make current human the infection source
    inline void specialR_0(){
        // At this point special TBV should have been given to everyone but the PEV to no-one
        effects[R_0_PEV].initialEfficacy = 1.0;
        effects[R_0_TBV].initialEfficacy = 0.0;
    }
#endif
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        effects & stream;
    }

private:
    /// Details for each deployed vaccine for this human
    typedef vector<PerEffectPerHumanVaccine> EffectList;
    EffectList effects;
};

}
}
#endif
