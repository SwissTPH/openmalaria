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

// Forward declare: doesn't need to be known about here
namespace scnXml {
class VaccineDescription;
}

namespace OM {
namespace Host {

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
    /// Set parameters from xml
    static void initParameters ();
    /// Free memory
    static void clearParameters ();

    /// True if any types of vaccine are in use.
    static bool anyVaccine;

private:
    /*! Common to all vaccine types. Number of vaccine doses that are given
     * either through EPI or as EPI Boosters. */
    static size_t _numberOfEpiDoses;

    /** Target age for EPI-like vaccination, in time steps. */
    static vector<int> targetAgeTStep;

    /// Preerythrocytic reduces h vaccine parameters
    static Vaccine PEV;
    /// Erythrocytic reduces y vaccine parameters
    static Vaccine BSV;
    /// Transmission blocking reduces k vaccine parameters
    static Vaccine TBV;

    //Non-static:
    Vaccine() : active(false), decay(1.0) {}

    /** Get the efficacy of the vaccine.
     *
     * @param numPrevDoses The number of prior vaccinations of the individual. */
    double getEfficacy (int numPrevDoses);

    /// True if this vaccine is in use
    bool active;

    /// exp(-Decay rate)
    double decay;
    /** Per-type initialization
     * @returns decay */
    void initVaccine (const scnXml::VaccineDescription* vd);

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

    /// Decay vaccines
    void update();

    /// Returns true if a continuous vaccine dose should be given.
    bool doCtsVaccination (int ageTSteps) {
        // Deployment is affected by previous missed doses and mass vaccinations,
        // unlike other continuous interventions; extra test:
        return _lastVaccineDose < (int)Vaccine::_numberOfEpiDoses
               && Vaccine::targetAgeTStep[_lastVaccineDose] == ageTSteps;
    }

    /** Update efficacies and the number of doses in this human. */
    void vaccinate();

    inline double getBSVEfficacy()const {
        return _BSVEfficacy;
    }
    inline double getPEVEfficacy()const {
        return _PEVEfficacy;
    }
    inline double getTBVEfficacy()const {
        return _TBVEfficacy;
    }

    /// Hack for R_0 deployment
    inline void removeTBV() {
        _TBVEfficacy = 0.0;
    }

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        _lastVaccineDose & stream;
        _BSVEfficacy & stream;
        _PEVEfficacy & stream;
        _TBVEfficacy & stream;
    }


private:
    /** Number of vaccine doses this individual has received.
     *
     * If an individual misses one EPI (continuous) vaccine dose, it's
     * intentional that they also miss following EPI doses (unless a timed mass
     * vaccination reintroduces them to the EPI schedule). */
    int _lastVaccineDose;
    //!Remaining efficacy of Blood-stage vaccines
    double _BSVEfficacy;
    //!Remaining efficacy of Pre-erythrocytic vaccines
    double _PEVEfficacy;
    //!Remaining efficacy of Transmission-blocking vaccines
    double _TBVEfficacy;
};

}
}
#endif
