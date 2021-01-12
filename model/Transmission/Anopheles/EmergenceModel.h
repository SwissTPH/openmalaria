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

#ifndef Hmod_Anopheles_EmergenceModel
#define Hmod_Anopheles_EmergenceModel

#include "Global.h"
#include "schema/interventions.h"
#include "util/SimpleDecayingValue.h"
#include "util/checkpoint_containers.h"
#include "util/vecDay.h"

#include <cmath>

#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {
namespace Anopheles {

using namespace std;
using namespace OM::util;

// forward declare to avoid circular dependency:
class MosqTransmission;

/** Part of vector anopheles model, giving emergence of adult mosquitoes from
 * water bodies.
 * 
 * This is an abstract class (i.e. interface). The following implementations
 * exist: FixedEmergence.
 */
class EmergenceModel
{
public:
    ///@brief Initialisation and destruction
    //@{
    EmergenceModel() {}
    virtual ~EmergenceModel() {}
    
    /** Latter part of AnophelesModel::init2.
     *
     * @param tsP_A P_A for initialisation; assumed constant when there are no interventions
     * @param tsP_df P_df for initialisation; assumed constant when there are no interventions
     * @param tsP_dff P_dff for initialisation; assumed constant when there are no interventions
     * @param EIRtoS_v multiplication factor to convert input EIR into required
     * @param transmission reference to MosqTransmission object
     * S_v. */
    virtual void init2(double tsP_dff, double initNv0FromSv, const vecDay<double>& forcedS_v, const vecDay<double>& mosqEmergeRate, const SimTime &mosqRestDuration) =0;
    
    /** Work out whether another interation is needed for initialisation and if
     * so, make necessary changes.
     *
     * @returns true if another iteration is needed. */
    virtual void initIterate (double factor, const vecDay<double>& mosqEmergeRate) = 0;
    //@}

    
    /** Model updates.
     * 
     * Returns the emergence for today, taking interventions like larviciding
     * into account, and updates some statistics (needed during
     * initialisation).
     * 
     * @param d0 Time of the start of the day-long update period
     * @param nOvipositing The number of adults which successfully
     * oviposited this/last time-step.
     * @param S_v Value of S_v for this day
     * @returns The number of adults emerging between the last simulated time
     * point and the one being calculated.
     */
    virtual double update(const SimTime &d0, const vecDay<double>& mosqEmergeRate, double nOvipositing) =0;

    virtual double getResAvailability() const =0;
    virtual double getResRequirements() const =0;

    template<class S>
    void operator& (S& stream) { checkpoint (stream); }
    
protected:
    // /** Return the proportion of emerging larvae to survive intervention
    //  * effects. Should be between 0 and 1. */
    // inline double interventionSurvival() const{ return emergenceSurvival; }
    
    virtual void checkpoint (istream& stream) =0;
    virtual void checkpoint (ostream& stream) =0;

private:
    /** Emergence rate of new mosquitoes, for every day of the year (N_v0).
     * 
     * Has annual periodicity: length is 365. First value (index 0) corresponds
     * to first day of year (1st Jan or something else if rebased). In 5-day
     * time-step model values at indecies 0 through 4 are used to calculate the
     * state at time-step 1.
     * 
     * Units: Animals per day.
     *
     * Should be checkpointed. */
    vecDay<double> mosqEmergeRate;
};

}
}
}
#endif
