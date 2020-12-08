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
    SimpleMPDEmergence(const scnXml::SimpleMPD& elt)
    {
        quinquennialOvipositing.assign (SimTime::fromYearsI(5), 0.0);
        invLarvalResources.assign (SimTime::oneYear(), 0.0);
        
        developmentDuration = SimTime::fromDays(elt.getDevelopmentDuration().getValue());
        if (!(developmentDuration > SimTime::zero()))
            throw util::xml_scenario_error("entomology.vector.simpleMPD.developmentDuration: "
                "must be positive");
        probPreadultSurvival = elt.getDevelopmentSurvival().getValue();
        if (!(0.0 <= probPreadultSurvival && probPreadultSurvival <= 1.0))
            throw util::xml_scenario_error("entomology.vector.simpleMPD.developmentSurvival: "
                "must be a probability (in range [0,1]");
        fEggsLaidByOviposit = elt.getFemaleEggsLaidByOviposit().getValue();
        if (!(fEggsLaidByOviposit > 0.0))
            throw util::xml_scenario_error("entomology.vector.simpleMPD.femaleEggsLaidByOviposit: "
                "must be positive");
        nOvipositingDelayed.assign (developmentDuration, 0.0);
    }

    /** Latter part of AnophelesModel::init2.
     *
     * @param tsP_A P_A for this time step.
     * @param tsP_df P_df for this time step.
     * @param tsP_dff P_dff for this time step.
     * @param EIRtoS_v multiplication factor to convert input EIR into required
     * @param transmission reference to MosqTransmission object
     * S_v. */
    void init2(double tsP_dff, double initNvFromSv, const vecDay<double>& forcedS_v, const vecDay<double>& mosqEmergeRate, const SimTime &mosqRestDuration)
    {
        // Initialise nOvipositingDelayed
        SimTime y1 = SimTime::oneYear();
        SimTime tau = mosqRestDuration;
        for( SimTime t = SimTime::zero(); t < developmentDuration; t += SimTime::oneDay() ){
            nOvipositingDelayed[mod_nn(t+tau, developmentDuration)] =
                tsP_dff * initNvFromSv * forcedS_v[t];
        }

        // Used when calculating invLarvalResources (but not a hard constraint):
        assert(tau+developmentDuration <= y1);
        for( SimTime t = SimTime::zero(); t < SimTime::oneYear(); t += SimTime::oneDay() ){
            double yt = fEggsLaidByOviposit * tsP_dff * initNvFromSv *
                forcedS_v[mod_nn(t + y1 - tau - developmentDuration, y1)];
            invLarvalResources[t] = (probPreadultSurvival * yt - mosqEmergeRate[t]) /
                (mosqEmergeRate[t] * yt);
        }
    }
    /** Work out whether another interation is needed for initialisation and if
     * so, make necessary changes.
     *
     * @returns true if another iteration is needed. */
    virtual void initIterate (double factor, const vecDay<double>& mosqEmergeRate)
    {
        vectors::scale (nOvipositingDelayed, factor);

        SimTime y1 = SimTime::oneYear(),
            y2 = SimTime::fromYearsI(2),
            y3 = SimTime::fromYearsI(3),
            y4 = SimTime::fromYearsI(4),
            y5 = SimTime::fromYearsI(5);
        assert(mosqEmergeRate.size() == y1);
        
        for( SimTime t = SimTime::zero(); t < y1; t += SimTime::oneDay() ){
            SimTime ttj = t - developmentDuration;
            // b · P_df · avg_N_v(t - θj - τ):
            double yt = fEggsLaidByOviposit * 0.2 * (
                quinquennialOvipositing[ttj + y1] +
                quinquennialOvipositing[ttj + y2] +
                quinquennialOvipositing[ttj + y3] +
                quinquennialOvipositing[ttj + y4] +
                quinquennialOvipositing[mod_nn(ttj + y5, y5)]);
            invLarvalResources[t] = (probPreadultSurvival * yt - mosqEmergeRate[t]) /
                (mosqEmergeRate[t] * yt);
        }
    }
    //@}
    
    virtual double update(const SimTime &d0, const vecDay<double>& mosqEmergeRate, double nOvipositing)
    {
        // Simple Mosquito Population Dynamics model: emergence depends on the
        // adult population, resources available, and larviciding.
        // See: A Simple Periodically-Forced Difference Equation Model for
        // Mosquito Population Dynamics, N. Chitnis, 2012. TODO: publish & link.
        SimTime d1 = d0 + SimTime::oneDay();

        double yt = fEggsLaidByOviposit * nOvipositingDelayed[mod_nn(d1, developmentDuration)];
        double emergence = probPreadultSurvival * yt / (1.0 + invLarvalResources[mod_nn(d0, SimTime::oneYear())] * yt);
        nOvipositingDelayed[mod_nn(d1, developmentDuration)] = nOvipositing;
        quinquennialOvipositing[mod_nn(d1, SimTime::fromYearsI(5))] = nOvipositing;
        return emergence;
    }
    ///@brief Interventions and reporting
    //@{
    double getResAvailability() const
    {
        //TODO: why offset by one time step? This is effectively getting the resources available on the last time step
        //TODO: only have to add one year because of offset
        SimTime start = sim::now() - SimTime::oneTS() + SimTime::oneYear();
        double total = 0;
        for( SimTime i = start, end = start + SimTime::oneTS(); i < end; i += SimTime::oneDay() ){
            SimTime dYear1 = mod_nn(i, SimTime::oneYear());
            total += 1.0 / invLarvalResources[dYear1];
        }
        return total / SimTime::oneTS().inDays();
    }

    double getResRequirements() const {
        return numeric_limits<double>::quiet_NaN();
    }
    //@}
    
protected:
    virtual void checkpoint (istream& stream){ (*this) & stream; }
    virtual void checkpoint (ostream& stream){ (*this) & stream; }
    
private:
    
    /// Checkpointing
    //Note: below comments about what does and doesn't need checkpointing are ignored here.
    template<class S>
    void operator& (S& stream) {
        quinquennialOvipositing & stream;
        developmentDuration & stream;
        probPreadultSurvival & stream;
        fEggsLaidByOviposit & stream;
        invLarvalResources & stream;
        nOvipositingDelayed & stream;
    }
    
    // -----  model parameters (loaded from XML)  -----
    
    /** Duration of development (time from egg laying to emergence) in days. */
    SimTime developmentDuration;
    
    /** Survival probability of a mosquito from egg to emergence in the absence
     * of density dependent mortality. */
    double probPreadultSurvival;
    
    /** Mean number of female eggs laid when a mosquito oviposites. */
    double fEggsLaidByOviposit;
    
    
    // -----  parameters (constant after initialisation)  -----
    
    /** As quinquennialS_v, but for N_v*P_df (units: animals). */
    vecDay<double> quinquennialOvipositing;
    
    
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
     * vecDay be checkpointed. */
    vecDay<double> invLarvalResources;
    
    /** Vector for storing values of nOvipositing for the last
     * developmentDuration time steps. Index 0 should correspond to
     * nOvipositing developmentDuration days before
     * get(0, dYear1, nOvipositing) is called. */
    vecDay<double> nOvipositingDelayed;
};

}
}
}
#endif
