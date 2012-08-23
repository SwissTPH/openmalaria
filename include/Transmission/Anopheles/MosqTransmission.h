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

#ifndef Hmod_Anopheles_Transmission
#define Hmod_Anopheles_Transmission

#include "Global.h"
#include "Monitoring/Survey.h"
#include "Transmission/Anopheles/EmergenceModel.h"
#include "schema/entomology.h"

#include <limits>
#include <boost/shared_ptr.hpp>

using boost::shared_ptr;

class MosqLifeCycleSuite;

namespace OM {
namespace Transmission {
namespace Anopheles {

// enumeration of gettable stats for cts out
enum VecStat { PA, PDF, PDIF, NV, OV, SV };

/** Encapsulates the central part of the Chitnis et al transmission model:
 * vector transmission of malaria.
 * 
 * This is only part of the model; the VectorModel class is largely just a
 * wrapper to support multiple mosquito species, and the AnophelesModel
 * class adds parameter initialisation and intervention support to this class
 * as well as translating between the (1- or 5-day) time-steps used by the
 * simulator and the 1-day time-step used by this model.
 */
class MosqTransmission {
public:
    ///@brief initialisation functions
    //@{
    MosqTransmission();
    
    /** Initialise parameters and variables.
     * 
     * This is only a fraction of parameter initialisation; see also
     * AnophelesModel::initialise. */
    void initialise ( const scnXml::AnophelesParams::LifeCycleOptional& lcOpt, const scnXml::AnophelesParams::SimpleMPDOptional& simpleMPDOpt, const scnXml::Mosq& mosq );
    
    /** (Re) allocate and initialise some state variables. Must be called
     * before model is run. */
    void initState ( double tsP_A, double tsP_df,
                     double initNvFromSv, double initOvFromSv,
                     const vector<double>& forcedS_v );
    
    /// Helper function for initialisation.
    void initIterateScale ( double factor );
    //@}
    
    /** Update by one day (may be called multiple times for 1 time-step update).
     * 
     * @param d The day whose state we are calculating.
     * @param tsP_A P_A for this time-step
     * @param tsP_df P_df for this time-step
     * @param tsP_dif P_dif for this time-step
     * @param printDebug Print some info to cerr
     * @returns S_v for the next time-step
     */
    double update( size_t d, double tsP_A, double tsP_df, double tsP_dif, bool isDynamic, bool printDebug );
    
    ///@brief Interventions and reporting
    //@{
    void uninfectVectors();
    //@}
    
    inline int getEIPDuration() const {
        return EIPDuration;
    }
    
    ///@brief Functions used in reporting
    //@{
    /// Reset per-time-step statistics before running time-step updates
    inline void resetTSStats() {
        timestep_N_v0 = 0.0;
    }
    /// Get mean emergence per day during last time-step
    inline double getLastN_v0 () const{
        return timestep_N_v0 / TimeStep::interval;
    }
    /// Get mean P_A/P_df/P_dif/N_v/O_v/S_v during last time-step
    /// @param vs PA, PDF, PDIF, NV, OV or SV
    double getLastVecStat ( VecStat vs ) const;
    
    inline double getMosqRestDuration() const {
        return mosqRestDuration;
    }
    
    inline double getResAvailability() const{
        return emergence->getResAvailability();
    }
    inline double getResRequirements() const{
        return emergence->getResRequirements();
    }
    
    /// Write some per-species summary information.
    void summarize (const string speciesName, Monitoring::Survey& survey) const;
    //@}
    
    
    /// Checkpointing
    //Note: below comments about what does and doesn't need checkpointing are ignored here.
    template<class S>
    void operator& (S& stream) {
        (*emergence) & stream;
        mosqRestDuration & stream;
        EIPDuration & stream;
        N_v_length & stream;
        P_A & stream;
        P_df & stream;
        P_dif & stream;
        N_v & stream;
        O_v & stream;
        S_v & stream;
        fArray & stream;
        ftauArray & stream;
        timestep_N_v0 & stream;
    }
    
    /** @brief Emergence model
     * 
     * Code to calculate emergence of mosquitoes from water bodies goes here. */
    shared_ptr<EmergenceModel> emergence;
    
private:
    // -----  parameters (constant after initialisation)  -----
    
    /** @brief Duration parameters for mosquito/parasite life-cycle
     * 
     * Currently these are all constant. In theory they could be made to vary
     * seasonally, based on a fixed periodic cycle, though some code and
     * possibly model changes would be needed to accomodate this.
     * 
     * All have units of days.
     *
     * Set in initialise function from XML data; no need to checkpoint. */
    //@{
    /** Duration of feeding cycle (equals duration of resting period) for
     * mosquito (τ).
     * Units: days. */
    int mosqRestDuration;

    /** Duration of the extrinsic incubation period (sporozoite development time)
    * (θ_s).
    * Units: Days.
    *
    * Doesn't need checkpointing. */
    int EIPDuration;
    
    /** N_v_length-1 is the number of previous days for which some parameters are
     * stored: P_A, P_df, P_dif, N_v, O_v and S_v. This is longer than some of
     * the arrays need to be, but simplifies code with no real impact.
     *
     * Should equal EIPDuration + mosqRestDuration to allow values up to
     * θ_s + τ - 1 days back, plus current day.
     *
     * Set by initialise; no need to checkpoint. */
    int N_v_length;
    //@}

    /// If less than this many mosquitoes remain infected, transmission is interrupted.
    double minInfectedThreshold;
    
    
    // -----  variable model state  -----
    
    /** @brief Variable arrays N_v_length long.
     *
     * P_A, P_df, P_dif, N_v, O_v and S_v are set in advancePeriod().
     *
     * Values at index ((d-1) mod N_v_length) are used to derive the state of
     * the population on day d. The state during days (t×(I-1)+1) through (t×I)
     * where t is TimeStep::simulation and I is TimeStep::interval is what
     * drives the transmission at time-step t.
     * 
     * These arrays should be checkpointed. */
    //@{
    /** Probability of a mosquito not finding a host one night. */
    vector<double> P_A;

    /** P_df and P_dif per-day.
     *
     * P_df is the probability of a mosquito finding a host and completing a
     * feeding cycle without being killed.
     *
     * P_dif is the probability of a mosquito finding a host, getting infected,
     * and successfully completing a feeding cycle.
     *
     * HOWEVER, if the initialisation phase is driven by an input EIR and not by
     * vector calculations, then during the initialisation phase, P_dif contains
     * the daily kappa values read from XML for validation purposes. */
    vector<double> P_df, P_dif;

    /** Numbers of host-seeking mosquitos each day
     * 
     * N_v is the total number of host-seeking mosquitoes, O_v is those seeking
     * and infected, and S_v is those seeking and infective (to humans). */
    vector<double> N_v, O_v, S_v;
    //@}

    ///@brief Working memory
    /** Used for calculations within advancePeriod. Only saved for optimisation.
     *
     * Used to calculate recursive functions f and f_τ in NDEMD eq 1.6, 1.7.
     * Values are recalculated each step; only fArray[0] and
     * ftauArray[0..mosqRestDuration] are stored across steps for optimisation
     * (reallocating each time they are needed would be slow).
     *
     * Length (fArray): EIPDuration - mosqRestDuration + 1 (θ_s - τ + 1)
     * Length (ftauArray): EIPDuration (θ_s)
     *
     * Don't need to be checkpointed, but some values need to be initialised. */
    //@{
    vector<double> fArray;
    vector<double> ftauArray;
    //@}
    
    /** Variables tracking data to be reported. */
    double timestep_N_v0;
    
    friend class ::MosqLifeCycleSuite;
};

}
}
}
#endif
