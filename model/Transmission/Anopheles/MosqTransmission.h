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

#ifndef Hmod_Anopheles_Transmission
#define Hmod_Anopheles_Transmission

#include "Global.h"
#include "Transmission/Anopheles/EmergenceModel.h"
#include "schema/entomology.h"

#include <limits>
#include <boost/shared_ptr.hpp>

class MosqLifeCycleSuite;

namespace OM {
namespace Transmission {
namespace Anopheles {
using util::vecDay2D;

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
    void initState ( double tsP_A, double tsP_df, double tsP_dff,
                     double initNvFromSv, double initOvFromSv,
                     const vecDay<double>& forcedS_v );
    
    /// Helper function for initialisation.
    void initIterateScale ( double factor );
    
    /** Set up the non-host-specific interventions. */
    inline void initVectorInterv( const scnXml::VectorSpeciesIntervention& elt, size_t instance ){
        emergence->initVectorInterv( elt, instance ); }
    //@}
    
    /** Update by one day (may be called multiple times for 1 time-step update).
     * 
     * @param d0 Time of the start of the day-long update period
     * @param tsP_A P_A for this time-step
     * @param tsP_df P_df for this time-step
     * @param tsP_dif P_dif for this time-step, per parasite genotype
     * @param tsP_dff P_dff for this time step
     * @param partialEIR Vector, per genotype; after calculation, the latest
     *  S_v values are multiplied by EIR_factor and added to this.
     * @param EIR_factor see parameter partialEIR
     */
    void update( SimTime d0, double tsP_A, double tsP_df,
                   const vector<double> tsP_dif, double tsP_dff,
                   bool isDynamic,
                   vector<double>& partialEIR, double EIR_factor );
    
    ///@brief Interventions and reporting
    //@{
    void uninfectVectors();
    //@}
    
    inline SimTime getEIPDuration() const {
        return EIPDuration;
    }
    
    ///@brief Functions used in reporting
    //@{
    /// Reset per-time-step statistics before running time-step updates
    inline void resetTSStats() {
        timeStep_N_v0 = 0.0;
    }
    /// Get mean emergence per day during last time-step
    inline double getLastN_v0 () const{
        return timeStep_N_v0 / SimTime::oneTS().inDays();
    }
    /// Get mean P_A/P_df/P_dif/N_v/O_v/S_v during last time-step
    /// @param vs PA, PDF, PDIF, NV, OV or SV
    double getLastVecStat( VecStat vs )const;
    
    inline SimTime getMosqRestDuration() const {
        return mosqRestDuration;
    }
    
    inline double getResAvailability() const{
        return emergence->getResAvailability();
    }
    inline double getResRequirements() const{
        return emergence->getResRequirements();
    }
    
    /// Write some per-species summary information.
    void summarize( size_t species )const;
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
        P_dff & stream;
        N_v & stream;
        O_v & stream;
        S_v & stream;
        //TODO: do we actually need to checkpoint these next three?
        fArray & stream;
        ftauArray & stream;
        uninfected_v & stream;
        timeStep_N_v0 & stream;
    }
    
    /** @brief Emergence model
     * 
     * Code to calculate emergence of mosquitoes from water bodies goes here. */
    boost::shared_ptr<EmergenceModel> emergence;
    
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
    SimTime mosqRestDuration;

    /** Duration of the extrinsic incubation period (sporozoite development time)
    * (θ_s).
    * Units: Days.
    *
    * Doesn't need checkpointing. */
    SimTime EIPDuration;
    
    /** N_v_length-1 is the number of previous days for which some parameters are
     * stored: P_A, P_df, P_dif, N_v, O_v and S_v. This is longer than some of
     * the arrays need to be, but simplifies code with no real impact.
     *
     * Should equal EIPDuration + mosqRestDuration to allow values up to
     * θ_s + τ - 1 days back, plus current day.
     *
     * Set by initialise; no need to checkpoint. */
    SimTime N_v_length;
    //@}
    
    /// If less than this many mosquitoes remain infected, transmission is interrupted.
    double minInfectedThreshold;
    
    
    // -----  variable model state  -----
    
    /** @brief Variable arrays N_v_length long.
     *
     * P_A, P_df, P_dif, N_v, O_v and S_v are set in advancePeriod() and have
     * values stored per day.
     * 
     * P_dif, O_v and S_v have a second index: the parasite genotype.
     *
     * Values at index ((d-1) mod N_v_length) are used to derive the state of
     * the population on day d. The state during days (t×I+1) through to ((t+1)×I)
     * where t is sim::ts0() and I is SimTime::oneTS().inDays() is what
     * drives the transmission at time-step t.
     * 
     * These arrays should be checkpointed. */
    //@{
    /** Probability of a mosquito not finding a host one night. */
    vecDay<double> P_A;
    
    /** P_df is the probability of a mosquito finding a host and completing a
     * feeding cycle without being killed. */
    vecDay<double> P_df;
    
    /** P_dif is the probability of a mosquito finding a host, getting
     * infected, and successfully completing a feeding cycle. */
    vecDay2D<double> P_dif;
    
    /** Like P_df but including fertility factors */
    vecDay<double> P_dff;
    
    /** Numbers of host-seeking mosquitos each day
     * 
     * N_v is the total number of host-seeking mosquitoes. */
    vecDay<double> N_v;
    
    /** Numbers of host-seeking mosquitos each day
     * 
     * O_v is the number of infected host-seeking mosquitoes, and S_v is the
     * number of infective (to humans) host-seeking mosquitoes. */
    vecDay2D<double> O_v, S_v;
    //@}

    ///@brief Working memory
    /** Used for calculations within advancePeriod. Only saved for optimisation.
     *
     * Used to calculate recursive functions f and f_τ in NDEMD eq 1.6, 1.7.
     * Values are recalculated each step; only fArray[0] and
     * ftauArray[0..mosqRestDuration] are stored across steps for optimisation
     * (reallocating each time they are needed would be slow).
     * 
     * uninfected_v[0] is not used.
     *
     * Length (fArray): EIPDuration - mosqRestDuration + 1 (θ_s - τ + 1)
     * Length (ftauArray): EIPDuration (θ_s)
     * Length (uninfected_v): N_v_length
     *
     * Don't need to be checkpointed, but some values need to be initialised. */
    //@{
    vecDay<double> fArray;
    vecDay<double> ftauArray;
    vecDay<double> uninfected_v;
    //@}
    
    /** Variables tracking data to be reported. */
    double timeStep_N_v0;
    
    friend class ::MosqLifeCycleSuite;
};

}
}
}
#endif
