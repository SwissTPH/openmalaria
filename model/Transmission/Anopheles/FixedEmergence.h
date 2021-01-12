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

#ifndef Hmod_Anopheles_FixedEmergence
#define Hmod_Anopheles_FixedEmergence

#include "Global.h"
#include "Transmission/Anopheles/EmergenceModel.h"
#include "schema/interventions.h"
#include <vector>
#include <limits>

#include "rotate.h"

namespace OM {
namespace Transmission {
namespace Anopheles {

using namespace std;
using namespace OM::util;

// forward declare to avoid circular dependency:
class MosqTransmission;

/** Part of vector anopheles model, giving emergence of adult mosquitoes from
 * water bodies. This model fits annual (periodic) sequence to produce the
 * desired EIR during warmup, then fixes this level of emergence for the rest
 * of the simulation.
 * 
 * Larviciding intervention directly scales the number of mosquitoes emerging
 * by a number, usually in the range [0,1] (but larger than 1 is also valid).
 */
class FixedEmergence : public EmergenceModel
{
public:
    /** Latter part of AnophelesModel::init2.
     *
     * @param tsP_A P_A for this time step.
     * @param tsP_df P_df for this time step.
     * @param tsP_dff P_dff for this time step.
     * @param EIRtoS_v multiplication factor to convert input EIR into required
     * @param transmission reference to MosqTransmission object
     * S_v. */
    virtual void init2(double tsP_dff, double initNv0FromSv, const vecDay<double>& forcedS_v, const vecDay<double>& mosqEmergeRate, const SimTime &mosqRestDuration){ }
    
    /** Work out whether another interation is needed for initialisation and if
     * so, make necessary changes.
     *
     * @returns true if another iteration is needed. */
    virtual void initIterate (double factor, const vecDay<double>& mosqEmergeRate) {}
    //@}
    
    virtual double update(const SimTime &d0, const vecDay<double>& mosqEmergeRate, double nOvipositing)
    {   
        // Get emergence at start of step:
        SimTime dYear1 = mod_nn(d0, SimTime::oneYear());
        // Simple model: fixed emergence scaled by larviciding
        return mosqEmergeRate[dYear1];
    }
    
    ///@brief Interventions and reporting
    //@{
    double getResAvailability() const {
        return numeric_limits<double>::quiet_NaN();
    }
    double getResRequirements() const {
        return numeric_limits<double>::quiet_NaN();
    }
    //@}
    
protected:
    virtual void checkpoint (istream& stream) {}
    virtual void checkpoint (ostream& stream) {}
};

}
}
}
#endif
