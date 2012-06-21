/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2012 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
 * TODO: most of this code is identical to that from the FixedEmergence model.
 * Use common base or even make this extend FixedEmergence?
 */
class SimpleMPDEmergence : public EmergenceModel
{
public:
    ///@brief Initialisation and destruction
    //@{
    /// Initialise and allocate memory
    SimpleMPDEmergence(const scnXml::SimpleMPD& elt);

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
    
    /// Return the emergence for today, taking interventions like larviciding
    /// into account.
    double get( size_t d, size_t dYear1, double nOvipositing );
    
    /// Store S_v for day d. Used by initIterate().
    void updateStats( size_t d, double tsP_dif, double S_v );
    
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
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
private:
    
    /// Checkpointing
    //Note: below comments about what does and doesn't need checkpointing are ignored here.
    template<class S>
    void operator& (S& stream) {
        EIRRotateAngle & stream;
        FSRotateAngle & stream;
        FSCoeffic & stream;
        mosqEmergeRate & stream;
        forcedS_v & stream;
        quinquennialS_v & stream;
        quinquennialOvipositing & stream;
        initNv0FromSv & stream;
        initNvFromSv & stream;
        initOvFromSv & stream;
        larvicidingEndStep & stream;
        larvicidingIneffectiveness & stream;
        developmentDuration & stream;
        probPreadultSurvival & stream;
        fEggsLaidByOviposit & stream;
        invLarvalResources & stream;
        nOvipositingDelayed & stream;
    }
    
    // -----  model parameters (loaded from XML)  -----
    
    /** Duration of development (time from egg laying to emergence) in days. */
    int developmentDuration;
    
    /** Survival probability of a mosquito from egg to emergence in the absence
     * of density dependent mortality. */
    double probPreadultSurvival;
    
    /** Mean number of female eggs laid when a mosquito oviposites. */
    double fEggsLaidByOviposit;
    
    
    // -----  parameters (constant after initialisation)  -----
    
    ///@brief Descriptions of transmission, used primarily during warmup
    //@{
    /** Summary of S_v over the last five years, used by vectorInitIterate to
     * calculate scaling factor.
     *
     * Length is 365 * 5. Checkpoint.
     *
     * Units: inoculations. */
    vector<double> quinquennialS_v;
    
    /** As quinquennialS_v, but for N_v*P_df (units: animals). */
    vector<double> quinquennialOvipositing;
    
    /** Conversion factor from forcedS_v to mosqEmergeRate.
     *
     *TODO: no longer true:
     * Also has another temporary use between initialise and setupNv0 calls:
     * "initOvFromSv" or  (ρ_O / ρ_S).
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
    vector<double> mosqEmergeRate;
    
    /** Resources for mosquito larvae (or rather 1 over resources); γ(t) in
     * model description.
     * 
     * Unlike model description, we allow special values 0 for no density
     * dependence and infinity for zero emergence.
     * 
     * Index t should correspond to the resources available to mosquitoes
     * emerging at t (i.e. same index in mosqEmergeRate).
     * 
     * Has annual periodicity: length is 365. First value (index 0) corresponds
     * to first day of year (1st Jan or something else if rebased). In 5-day
     * time-step model values at indecies 0 through 4 are used to calculate the
     * state at time-step 1.
     * 
     * Units: 1 / animals per day.
     *
     * Should be checkpointed. */
    vector<double> invLarvalResources;
    
    /** Vector for storing values of nOvipositing for the last
     * developmentDuration time steps. Index 0 should correspond to
     * nOvipositing developmentDuration days before
     * get(0, dYear1, nOvipositing) is called. */
    vector<double> nOvipositingDelayed;
};

}
}
}
#endif
