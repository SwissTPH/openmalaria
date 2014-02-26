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
#ifndef Hmod_Vaccine
#define Hmod_Vaccine

#include "schema/interventions.h"
#include "interventions/Interfaces.hpp"
#include "Global.h"
#include "util/DecayFunction.h"
#include <boost/shared_ptr.hpp>

namespace OM {
namespace Host {
    class Human;
}
namespace interventions {
    using util::DecayFunction;
    using util::DecayFuncHet;
    using boost::shared_ptr;
    
/** Vaccine intervention parameters.
 *
 * Used to represent PEV, BSV and TBV vaccines.
 * Each that has a descriptor is applied
 * simultaneously by a continuous or timed intervention strategy (no
 * way to determine which are used).
 *
 * All parameters (inc. non-static) are only set by initParameters(). */
class Vaccine {
public:
    enum Types { PEV, BSV, TBV, NumVaccineTypes };

    Vaccine(const scnXml::VaccineDescription& vd, Types type);
    
private:
    /** Get the initial efficacy of the vaccine.
     *
     * @param numPrevDoses The number of prior vaccinations of the individual. */
    double getInitialEfficacy (size_t numPrevDoses) const;

    inline static const Vaccine& getParams( Types type ){
        assert( params[type] != 0 );
        return *params[type];
    }
    
    /** @brief Vaccine static parameters
     * 
     * Each instance is either null or points to data for the vaccine type.
     *
     * No memory management (only leak is at exit which OS deals with). */
    static Vaccine* params[NumVaccineTypes];
    
    /** Until the monitoring system is updated, only one type of vaccination
     * should be reported. Which to report is described here. */
    static Types reportType;

    /// Function representing decay of effect
    shared_ptr<DecayFunction> decayFunc;

    /* Vaccine type specific parameters
     * Initial mean efficacy, definition depends on vaccine type */
    vector<double> initialMeanEfficacy;
    // Distribution of efficacies among individuals, parameter to sample from beta dist.
    double efficacyB;

    friend class PerHumanVaccine;
    friend class PerEffectPerHumanVaccine;
};

/** Per vaccine effect (type), per human details. */
class PerEffectPerHumanVaccine {
public:
    //Note: this constructor is only for checkpointing
    PerEffectPerHumanVaccine();
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        numDosesAdministered & stream;
        timeLastDeployment & stream;
        initialEfficacy & stream;
        hetSample & stream;
    }
    
private:
    PerEffectPerHumanVaccine( const Vaccine& params );
    
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
    DecayFuncHet hetSample;
    
    friend class PerHumanVaccine;
};

/** Per-human vaccine code. */
class PerHumanVaccine {
public:
    PerHumanVaccine() {
        types[0] = types[1] = types[2] = 0;
    }
    ~PerHumanVaccine(){
        for( size_t i = 0; i < Vaccine::NumVaccineTypes; ++i )
            if( types[i] != 0 )
                delete types[i];
    }
    
    /** Get one minus the efficacy of the vaccine (1 for no effect, 0 for full effect). */
    double getFactor( Vaccine::Types type )const;
    
    void possiblyVaccinate( const Host::Human& human,
                           Deployment::Method method,
                           Vaccine::Types type,
                           interventions::VaccineLimits vaccLimits );

#if 0
    /// Hack for R_0 experiment: make current human the infection source
    inline void specialR_0(){
        assert( Vaccine::types[Vaccine::PEV] != 0 && Vaccine::types[Vaccine::TBV] != 0 );
        types[Vaccine::PEV].initialEfficacy = 1.0;
        types[Vaccine::TBV].initialEfficacy = 0.0;
    }
#endif
    
    /// Checkpointing
    void operator& (ostream& stream) {
        for( size_t i = 0; i < Vaccine::NumVaccineTypes; ++i ){
            if( types[i] == 0 ){
                false & stream;
            }else{
                true & stream;
                (*types[i]) & stream;
            }
        }
    }
    void operator& (istream& stream) {
        for( size_t i = 0; i < Vaccine::NumVaccineTypes; ++i ){
            bool hasEffect;
            hasEffect & stream;
            if( hasEffect ){
                types[i] = new PerEffectPerHumanVaccine();
                (*types[i]) & stream;
            }else{
                types[i] = 0;
            }
        }
    }

private:
    /// Details for each vaccine type
    PerEffectPerHumanVaccine* types[Vaccine::NumVaccineTypes];
};

}
}
#endif
