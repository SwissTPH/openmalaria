/*
 This file is part of OpenMalaria.

 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

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
#ifndef Hmod_Vaccine
#define Hmod_Vaccine

#include "Global.h"
#include "util/DecayFunction.h"
#include "schema/interventions.h"
#include <boost/shared_ptr.hpp>

namespace OM {
namespace Host {
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
    // Static:
    /// Set parameters from xml (only called if vaccines are used)
    static void init (const scnXml::Vaccine& xmlVaccine);
    
private:
    /*! Common to all vaccine types. Number of vaccine doses that are given
     * either through EPI or as EPI Boosters. */
    static size_t _numberOfEpiDoses;

    /** Target age for EPI-like vaccination, in time steps. */
    static vector<TimeStep> targetAgeTStep;

    /// Preerythrocytic reduces h vaccine parameters
    static Vaccine PEV;
    /// Erythrocytic reduces y vaccine parameters
    static Vaccine BSV;
    /// Transmission blocking reduces k vaccine parameters
    static Vaccine TBV;

    //Non-static:
    Vaccine() : active(false), decayFunc(DecayFunction::makeConstantObject()), efficacyB(1.0) {}
    
    /** Per-type initialization
     * @returns decay */
    void initVaccine (const scnXml::VaccineDescription* vd);

    /** Get the efficacy of the vaccine.
     *
     * @param numPrevDoses The number of prior vaccinations of the individual. */
    double getEfficacy (int numPrevDoses);

    /// True if this vaccine is in use
    bool active;

    /// Function representing decay of effect
    shared_ptr<DecayFunction> decayFunc;

    /* Vaccine type specific parameters
     * Initial mean efficacy, definition depends on vaccine type */
    vector<double> initialMeanEfficacy;
    // Distribution of efficacies among individuals, parameter to sample from beta dist.
    double efficacyB;

    friend class PerHumanVaccine;
};

/** Per-human vaccine code. */
class PerHumanVaccine {
public:
    PerHumanVaccine();

    /// Returns true if a continuous vaccine dose should be given.
    bool doCtsVaccination (TimeStep ageTSteps) {
        // Deployment is affected by previous missed doses and mass vaccinations,
        // unlike other continuous interventions; extra test:
        return _lastVaccineDose < (int)Vaccine::_numberOfEpiDoses
               && Vaccine::targetAgeTStep[_lastVaccineDose] == ageTSteps;
    }

    /** Update efficacies and the number of doses in this human. */
    void vaccinate(TimeStep now);
    /// Has been vaccinated within considered effective duration?
    bool hasProtection(TimeStep maxInterventionAge)const;

    inline double getPEVEfficacy()const {
        return _initialPEVEfficacy * Vaccine::PEV.decayFunc->eval( TimeStep::simulation1() - _timeLastVaccine, hetSamplePEV );
    }
    inline double getBSVEfficacy()const {
        return _initialBSVEfficacy * Vaccine::BSV.decayFunc->eval( TimeStep::simulation1() - _timeLastVaccine, hetSampleBSV );
    }
    inline double getTBVEfficacy()const {
        return _initialTBVEfficacy * Vaccine::TBV.decayFunc->eval( TimeStep::simulation1() - _timeLastVaccine, hetSampleTBV );
    }
    
    /// @brief Hacks for R_0 deployment
    //@{
    inline void setInitialPEV( double effic ) {
	_initialPEVEfficacy = effic;
    }
    inline void setInitialTBV( double effic ) {
        _initialTBVEfficacy = effic;
    }
    //@}

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        _lastVaccineDose & stream;
        _timeLastVaccine & stream;
        _initialBSVEfficacy & stream;
        _initialPEVEfficacy & stream;
        _initialTBVEfficacy & stream;
        hetSamplePEV & stream;
        hetSampleBSV & stream;
        hetSampleTBV & stream;
    }


private:
    /** Number of vaccine doses this individual has received.
     *
     * If an individual misses one EPI (continuous) vaccine dose, it's
     * intentional that they also miss following EPI doses (unless a timed mass
     * vaccination reintroduces them to the EPI schedule). */
    int _lastVaccineDose;
    /// Timestep of last vaccination
    TimeStep _timeLastVaccine;
    //!Remaining efficacy of Pre-erythrocytic vaccines
    double _initialPEVEfficacy;
    //!Remaining efficacy of Blood-stage vaccines
    double _initialBSVEfficacy;
    //!Remaining efficacy of Transmission-blocking vaccines
    double _initialTBVEfficacy;
    
    DecayFuncHet hetSamplePEV;
    DecayFuncHet hetSampleBSV;
    DecayFuncHet hetSampleTBV;
};

}
}
#endif
