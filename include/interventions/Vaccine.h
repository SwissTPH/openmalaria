/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2013 Swiss Tropical and Public Health Institute 
 * Copyright (C) 2005-2013 Liverpool School Of Tropical Medicine
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
#include "interventions/Interventions.h"
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

    // ———  static  ———
    /// Special for R_0: check is set up correctly or throw xml_scenario_error
    static void verifyEnabledForR_0 ();
    
    /** Per-type initialization */
    inline static void init( const scnXml::VaccineDescription& vd, Types type ){
        types[type].initVaccine( vd, type );
    }
    
    /// Set schedule. Needed for correct EPI deployment.
    /// TODO(vaccines): a model of how vaccine booster shots work would allow this to be
    /// moved to intervention deployment.
    inline static void initSchedule( Types type, const scnXml::ContinuousList::DeploySequence& schedule ){
        types[type].initScheduleInternal( schedule );
    }
    
private:
    // ———  non-static  ———
    
    Vaccine() : active(false), decayFunc(DecayFunction::makeConstantObject()), efficacyB(1.0) {}
    
    void initVaccine (const scnXml::VaccineDescription& vd, Types type);
    
    void initScheduleInternal( const scnXml::ContinuousList::DeploySequence& schedule );
    
    /** Get the initial efficacy of the vaccine.
     *
     * @param numPrevDoses The number of prior vaccinations of the individual. */
    double getInitialEfficacy (size_t numPrevDoses);

    /** @brief Three types of vaccine.
     * 
     * TODO(vaccines): multiple descriptions should be allowed for each type. */
    static Vaccine types[NumVaccineTypes];
    
    /** Until the monitoring system is updated, only one type of vaccination
     * should be reported. Which to report is described here. */
    static Types reportType;

    /// True if this vaccine is in use
    bool active;

    /** Target age for EPI-like vaccination, in time steps. */
    vector<TimeStep> targetAgeTStep;

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
    PerEffectPerHumanVaccine();
    
    double getEfficacy( Vaccine::Types type ) const;
    
    /// Returns true if this individual should get a vaccine dose via EPI
    bool getsEPIVaccination( Vaccine::Types type, TimeStep ageTSteps ) const;

    /** Update efficacies and the number of doses in this human. */
    void vaccinate( Deployment::Method method,
                    Vaccine::Types type );
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        numDosesAdministered & stream;
        timeLastDeployment & stream;
        initialEfficacy & stream;
        hetSample & stream;
    }
    
private:
    PerEffectPerHumanVaccine(Vaccine::Types type );
    
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
    PerHumanVaccine();
    
    inline double getEfficacy( Vaccine::Types type )const{
        return types[type].getEfficacy( type );
    }
    
    void vaccinate( const Host::Human& human,
                           Deployment::Method method,
                           Vaccine::Types type );

    /// Hack for R_0 experiment: make current human the infection source
    inline void specialR_0(){
        assert( Vaccine::types[Vaccine::PEV].active && Vaccine::types[Vaccine::TBV].active );
        types[Vaccine::PEV].initialEfficacy = 1.0;
        types[Vaccine::TBV].initialEfficacy = 0.0;
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        for( size_t i = 0; i < Vaccine::NumVaccineTypes; ++i ){
            types[i] & stream;
        }
    }

private:
    /// Details for each vaccine type
    PerEffectPerHumanVaccine types[Vaccine::NumVaccineTypes];
};

}
}
#endif
