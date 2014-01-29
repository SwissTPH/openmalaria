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

#ifndef Hmod_WithinHost_Immunity
#define Hmod_WithinHost_Immunity

#include "Global.h"
#include "WithinHost/WHInterface.h"

#include <list>

using namespace std;

class UnittestUtil;

namespace OM {
namespace WithinHost {

/**
 * Immunity code and base class for all current P. falciparum models.
 */
class WHFalciparum : public WHInterface {
public:
    /// @brief Static methods
    //@{
    /// Initialise static parameters
    static void init();
    //@}

    /// @brief Constructors, destructors and checkpointing functions
    //@{
    WHFalciparum();
    virtual ~WHFalciparum();

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }
    //@}

    /// Create a new infection within this human
    virtual void importInfection() =0;
    /** Conditionally clears all infections. Not used with the PK/PD model.
     *
     * If IPT isn't present, it just calls clearAllInfections(); otherwise it
     * uses IPT code to determine whether to clear all infections or do nothing
     * (isSevere is only used in the IPT case). */
    virtual void clearInfections (bool isSevere);

    /** Medicate drugs (wraps drug's medicate).
     *
     * @param drugAbbrev	abbrevation of drug name (e.g. CQ, MF)
     * @param qty		Quantity of drug to administer in mg
     * @param time		Time relative to beginning of timestep to medicate at, in days (less than 1 day)
     * @param duration Duration in days. 0 or NaN indicate oral treatment.
     * @param bodyMass	Weight of human in kg
     */
    virtual void medicate(string drugAbbrev, double qty, double time, double duration, double bodyMass) {}


    inline double getCumulativeh() const {
        return _cumulativeh;
    }
    inline double getCumulativeY() const {
        return _cumulativeY;
    }

protected:
    ///@brief Immunity model
    //@{
    virtual void immunityPenalisation();

    /** Updates for the immunity model − assumes _cumulativeh and _cumulativeY
     * have already been incremented.
     *
     * Applies decay of immunity against asexual blood stages, if present. */
    void updateImmuneStatus();

    //!innate ability to control parasite densities
    double _innateImmSurvFact;

    /** Number of infections received since birth. */
    double _cumulativeh;
    //!Cumulative parasite density since birth
    double _cumulativeY;
    //!cumulativeY from previous timestep
    double _cumulativeYlag;
    //@}

    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

    ///@brief Static parameters, set by init()
    //@{
//Standard dev innate immunity for densities
    static double sigma_i;
// contribution of parasite densities to acquired immunity in the presence of fever
    static double immPenalty_22;
    /*
      Remaining immunity against asexual parasites(after time step, each of 2 components y and h)
      This variable decays the effectors cumulativeH and cumulativeY in a way that their
      effects on densities (1-Dh and 1-Dy) decay exponentially.
    */
    static double asexImmRemain;
    /*
      Remaining immunity against asexual parasites(after each time step, each of 2 components y and h)
      This variable decays the effectors cumulativeH and cumulativeY exponentially.
    */
    static double immEffectorRemain;
    //@}

    friend class ::UnittestUtil;
};

}
}
#endif
