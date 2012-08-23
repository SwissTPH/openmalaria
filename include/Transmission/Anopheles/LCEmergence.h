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

#ifndef Hmod_Anopheles_LCEmergence
#define Hmod_Anopheles_LCEmergence

#include "Global.h"
#include "Transmission/Anopheles/EmergenceModel.h"
//TODO: should LifeCycle be integrated here?
#include "Transmission/Anopheles/LifeCycle.h"
#include "schema/interventions.h"
#include <vector>

namespace OM {
namespace Transmission {
namespace Anopheles {

using namespace std;

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
class LCEmergence : public EmergenceModel
{
public:
    ///@brief Initialisation and destruction
    //@{
    /// Initialise and allocate memory
    LCEmergence();
    
    /** Called to initialise life-cycle parameters from XML data. */
    void initLifeCycle(const scnXml::LifeCycle& lcData);
    
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
    double getResAvailability() const{
        return lcParams.getResAvailability();
    }
    double getResRequirements() const{
        return lifeCycle.getResRequirements( lcParams );
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
        forcedS_v & stream;
        quinquennialP_dif & stream;
        initNv0FromSv & stream;
        initNvFromSv & stream;
        initOvFromSv & stream;
        initialP_A & stream;
        initialP_df & stream;
        lcParams & stream;
        lifeCycle & stream;
        larvicidingEndStep & stream;
        larvicidingIneffectiveness & stream;
    }
    
    // -----  parameters (constant after initialisation)  -----
    
    ///@brief Descriptions of transmission, used primarily during warmup
    //@{
    /** Summary of P_dif over the last five years, used by vectorInitIterate to
     * estimate larvalResources.
     *
     * Length is 365 * 5. Checkpoint.
     * 
     * NOTE: technically, only a fifth as many values need to be stored since
     * this only changes every five days. But that makes life more complicated.
     */
    vector<double> quinquennialP_dif;
    //@}
    
    ///@brief More stuff (init only?)
    //@{
    /** Conversion factor from forcedS_v to mosqEmergeRate.
     *
     *TODO: no longer true:
     * Also has another temporary use between initialise and setupNv0 calls:
     * "initOvFromSv" or  (ρ_O / ρ_S).
     *
     * Should be checkpointed. */
    double initNv0FromSv;       ///< ditto
    
    /** Values of P_A and P_df from initial population age structure. In theory
     * these values are constant until interventions start to affect mosquitoes
     * unless age structure varies due to low pop size or very high death
     * rates. */
    double initialP_A, initialP_df;
    //@}
    
    /// Parameters for life-cycle (excluding parasite transmission)
    LifeCycleParams lcParams;
    
    /// Mosquito life-cycle state
    LifeCycle lifeCycle;
    
    friend class ::MosqLifeCycleSuite;
};

}
}
}
#endif
