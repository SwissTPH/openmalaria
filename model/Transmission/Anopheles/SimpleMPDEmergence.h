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

#ifndef Hmod_Anopheles_SimpleMPDEmergence
#define Hmod_Anopheles_SimpleMPDEmergence

#include "Global.h"
#include "Transmission/Anopheles/EmergenceModel.h"
#include "schema/interventions.h"
#include <vector>
#include <limits>

namespace OM {
namespace Transmission {
namespace Anopheles {

using namespace std;

// forward declare to avoid circular dependency:
class MosqTransmission;

/** Part of vector anopheles model, giving emergence of adult mosquitoes from
 * water bodies. This model fits annual (periodic) sequence to produce the
 * desired EIR during warmup, then calculates larval resources (space) needed
 * to reproduce this emergence according to a simple model.
 * 
 * Larviciding intervention directly scales the number of mosquitoes emerging
 * by a number, usually in the range [0,1] (but larger than 1 is also valid).
 * Simple mosquito population dynamics model ensures reduction in adult
 * mosquito numbers affects emergence.
 * 
 * TODO(vec lifecycle): most of this code is identical to that from the FixedEmergence model.
 * Use common base or even make this extend FixedEmergence?
 */
class SimpleMPDEmergence : public EmergenceModel
{
public:
    ///@brief Initialisation and destruction
    //@{
    /// Initialise and allocate memory
    /// 
    /// @param species Index of species in list (numerical identifier, from 0 to num-species - 1)
    SimpleMPDEmergence(const scnXml::SimpleMPD& elt, size_t species);
    
    /// Static function which needs to be called exactly once before init2 or initIterate
    /// is called for each species (but after the SimpleMPDEmergence constructor).
    static void initShared();

    /** Latter part of AnophelesModel::init2.
     *
     * @param tsP_A P_A for this time step.
     * @param tsP_df P_df for this time step.
     * @param EIRtoS_v multiplication factor to convert input EIR into required
     * @param transmission reference to MosqTransmission object
     * S_v. */
    void init2( double tsP_A, double tsP_df, double EIRtoS_v, MosqTransmission& transmission );
    
    /** Work out whether another interation is needed for initialisation and if
     * so, make necessary changes.
     *
     * @returns true if another iteration is needed. */
    bool initIterate (MosqTransmission& transmission);
    //@}
    
    virtual double update( SimTime d0, double nOvipositing, double S_v );
    
    ///@brief Interventions and reporting
    //@{
    double getResAvailability() const;
    double getResRequirements() const {
        return numeric_limits<double>::quiet_NaN();
    }
    //@}
    
protected:
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
private:
    /// Checkpointing
    //Note: below comments about what does and doesn't need checkpointing are ignored here.
    template<class S>
    void operator& (S& stream) {
        staticCheckpoint(species, resType, stream);
        mosqEmergeRate & stream;
        quinquennialS_v & stream;
        quinquennialOvipositing & stream;
        initNv0FromSv & stream;
        probPreadultSurvival & stream;
        species & stream;
        resType & stream;
    }
    
    static void staticCheckpoint(size_t species, size_t resType, ostream& stream);
    static void staticCheckpoint(size_t species, size_t resType, istream& stream);
    
    // -----  model parameters (loaded from XML)  -----
    
    /** Survival probability of a mosquito from egg to emergence in the absence
     * of density dependent mortality. */
    double probPreadultSurvival;
    
    
    // -----  parameters (constant after initialisation)  -----
    
    /// Index of this mosquito species (see VectorModel).
    size_t species;
    
    /// Type of larval resources; used as key in internal LR namespace
    size_t resType;
    
    ///@brief Descriptions of transmission, used primarily during warmup
    //@{
    /** Summary of S_v over the last five years, used by vectorInitIterate to
     * calculate scaling factor.
     *
     * Length is 365 * 5. Checkpoint.
     *
     * Units: inoculations. */
    vecDay<double> quinquennialS_v;
    
    /** As quinquennialS_v, but for N_v*P_df (units: animals). */
    vecDay<double> quinquennialOvipositing;
    
    /** Conversion factor from forcedS_v to mosqEmergeRate.
     *
     * Should be checkpointed. */
    double initNv0FromSv;       ///< ditto
    //@}
    
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
