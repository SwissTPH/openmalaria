/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef Hmod_SpeciesModel
#define Hmod_SpeciesModel

#include "Global.h"
#include "Monitoring/Survey.h"
#include "Transmission/Vector/PerHost.h"
#include "Transmission/PerHost.h"
#include "Transmission/Vector/MosquitoLifeCycle.h"
#include <list>
#include <vector>

namespace OM {
namespace Host {
class Human;
}
namespace Transmission {
namespace Vector {

using namespace std;

/** Per-species part for vector transmission model.
 *
 * Data in this class is specific to the "species" of anopheles mosquito, where
 * species is used in a relaxed way to mean any variation of anopheles
 * mosquito, not just those types formally recognised as distinct species.
 *
 * A list of this class type is used by the VectorModel class to hold
 * (potentially) species-specific per-population data.
 *
 * Instead of storing static variables in this class, store them in the
 * VectorModel class.
 *
 * Variable names largely come from Nakul Chitnis's paper:
 * "A mathematical model for the dynamics of malaria in mosquitoes feeding on
 * a heterogeneous host population" (3rd Oct. 2007). */
class SpeciesModel
{
public:
    ///@brief Initialisation and destruction
    //@{
    SpeciesModel (const ITNParams* baseITNParams) :
            humanBase(baseITNParams),
            partialEIR(0.0),
            larvicidingEndStep (TimeStep::future),
            larvicidingIneffectiveness (1.0),
            timestep_N_v0(0.0)
    {}

    /** Called to initialise variables instead of a constructor. At this point,
     * the size of the human population is known but that population has not
     * yet been constructed. Called whether data is loaded from a check-point
     * or not.
     *
     * @param anoph Data structure from XML to use
     * @param initialisationEIR In/out parameter: TransmissionModel::initialisationEIR
     * @param nonHumanHostPopulations Size of each non-human population
     * @param populationSize Size of human population (assumed constant)
     */
    string initialise (const scnXml::AnophelesParams& anoph,
                       vector<double>& initialisationEIR,
                       map<string, double>& nonHumanHostPopulations,
                       int populationSize);

    /** Scale the internal EIR representation by factor; used as part of
     * initialisation. */
    void scaleEIR( double factor );

    /** Initialisation which must wait until a human population is available.
     * This is only called when a checkpoint is not loaded.
     *
     * @param sIndex Index in VectorModel::species of this class.
     * @param population The human population
     * @param populationSize Number of humans (use instead of population.size())
     *
     * Can only usefully run its calculations when not checkpointing, due to
     * population not being the same when loaded from a checkpoint. */
    void setupNv0 (size_t sIndex,
                   const std::list<Host::Human>& population,
                   int populationSize,
                   double invMeanPopAvail);

    /** Return base-line human parameters for the mosquito. */
    inline const Vector::PerHostBase& getHumanBaseParams () {
        return humanBase;
    }

    /** Set up intervention descriptions for humans, for this anopheles species. */
    inline void setITNDescription (const ITNParams& params,
                        const scnXml::ITNDescription::AnophelesParamsType& elt,
                        double proportionUse) {
        humanBase.setITNDescription (params, elt, proportionUse);
    }
    /** Set up intervention descriptions for humans, for this anopheles species. */
    inline void setIRSDescription (const scnXml::IRSDescription& irsDesc) {
        humanBase.setIRSDescription (irsDesc);
    }
    /** Set up intervention descriptions for humans, for this anopheles species. */
    inline void setVADescription (const scnXml::BaseInterventionDescription& vaDesc) {
        humanBase.setVADescription (vaDesc);
    }

    /** Work out whether another interation is needed for initialisation and if
     * so, make necessary changes.
     *
     * @returns true if another iteration is needed. */
    bool vectorInitIterate ();
    //@}
    
    
    ///@brief Functions called as part of usual per-timestep operations
    //@{
    /** Called per time-step. Does most of calculation of EIR.
     *
     * @param population The human population; so we can sum up availability and
     *  infectiousness.
     * @param populationSize Number of humans
     * @param sIndex Index of the type of mosquito in per-type/species lists.
     * @param isDynamic True to use full model; false to drive model from current contents of S_v.
     * @param invMeanPopAvail 1 over mean population availability relative to an adult.
     */
    void advancePeriod (const std::list<Host::Human>& population,
                        int populationSize,
                        size_t sIndex,
                        bool isDynamic,
                        double invMeanPopAvail);

    /** Returns the EIR calculated by advancePeriod().
     *
     * Could be extended to allow input EIR driven initialisation on a per-species
     * level instead of the whole simulation, but that doesn't appear worth doing.
     *
     * @param sIndex Index of this in VectorModel::species
     * @param host PerHost of the human requesting this EIR. */
    double calculateEIR (size_t sIndex, Transmission::PerHost& host) {
        if ( partialEIR != partialEIR ) {
            cerr<<"partialEIR is not a number; "<<sIndex<<endl;
        }
        /* Calculates EIR per individual (hence N_i == 1).
         *
         * See comment in SpeciesModel.advancePeriod for method. */
        return partialEIR
               * host.entoAvailabilityHetVecItv (humanBase, sIndex)
               * host.probMosqBiting(humanBase, sIndex);        // probability of biting, once commited
    }
    //@}


    ///@brief Functions called to deploy interventions
    //@{
    void intervLarviciding (/*const scnXml::LarvicidingAnopheles&*/);

    void uninfectVectors();
    //@}

    
    ///@brief Functions used in reporting
    //@{
    /// Get emergence during last time-step
    inline double getLastN_v0 () const{
        return timestep_N_v0;
    }
    enum VecStat { PA, PDF, PDIF, NV, OV, SV };
    /// Get P_A/P_df/P_dif/N_v/O_v/S_v during last time-step
    /// @param vs PA, PDF, PDIF, NV, OV or SV
    double getLastVecStat ( VecStat vs ) const;
    inline double getResAvailability() const{
        return lcParams.getResAvailability();
    }
    inline double getResRequirements() const{
        return lcModel.getResRequirements( lcParams );
    }

    /// Write some per-species summary information.
    void summarize (const string speciesName, Monitoring::Survey& survey) const;
    //@}
    

    /// Checkpointing
    //Note: below comments about what does and doesn't need checkpointing are ignored here.
    template<class S>
    void operator& (S& stream) {
        mosqSeekingDeathRate & stream;
        mosqSeekingDuration & stream;
        mosqRestDuration & stream;
        EIPDuration & stream;
        probMosqSurvivalOvipositing & stream;
        EIRRotateAngle & stream;
        FSRotateAngle & stream;
        FSCoeffic & stream;
#if 0
        mosqEmergeRate & stream;
#endif
        forcedS_v & stream;
        quinquennialS_v & stream;
        initNv0FromSv & stream;
        initNvFromSv & stream;
        N_v_length & stream;
        P_A & stream;
        P_df & stream;
        P_dif & stream;
        N_v & stream;
        O_v & stream;
        S_v & stream;
        fArray & stream;
        ftauArray & stream;
        partialEIR & stream;
        larvicidingEndStep & stream;
        larvicidingIneffectiveness & stream;
        timestep_N_v0 & stream;
    }


private:
    ///@brief Initialisation helper functions
    //@{
    /** Calculate availability rate of hosts (α_i) and death rate while seeking
     * (µ_vA)
     *
     * Documentation: "Parameter Values for Transmission model"
     * (Chitnis, Smith and Schapira, 4.3.2010)
     * 
     * @param anoph Data from XML
     * @param nonHumanHostPopulations Size of each non-human population
     * @param populationSize Size of the human population (assumed constant)
     */
    void initAvailability(
        const scnXml::AnophelesParams& anoph,
        map<string, double>& nonHumanHostPopulations,
        int populationSize);

    /** Calculates the human ento availability
     * 
     * Reference: Parameter Values for Transmission Model, Chitnis et al,
     * September 2010 eqn (26).
     * 
     * @param N_i Human/non-human population size
     * @param P_A Probability of mosquito not dying or finding a host while
     *  seeking on a given night
     * @param P_Ai Probability of mosquito finding a human/non-human host of
     *  type i while seeking on a given night
     * @return α_i, the rate at which mosquitoes encounter hosts of type i
     *  while seeking
     */
    double calcEntoAvailability(double N_i, double P_A, double P_Ai);
    
    /** Called by initialise function to init variables directly related to EIR
     * 
     * @param anoph Data from XML
     * @param initialisationEIR In/out parameter: TransmissionModel::initialisationEIR
     */
    void initEIR(
        const scnXml::AnophelesParams& anoph,
        vector<double>& initialisationEIR);
    //@}
    
    
    // -----  Variable/constant parameters  -----
    
    /** Baseline parameters which may be varied per human host. The primary
     * reason for wrapping these parameters in a struct is that these are the
     * parameters which need to be passed to the PerHost code
     * during calculations.
     *
     * Includes model parameters which may be varied per-individual to account
     * for interventions and innate resistances, and intervention effect
     * descriptions.
     *
     * Read from XML by initialise; no need to checkpoint. */
    Vector::PerHostBase humanBase;
    
    
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

    /** Duration of host-seeking per day; the maximum fraction of a day that a
     * mosquito would spend seeking (θ_d). */
    double mosqSeekingDuration;
    //@}
    
    
    /// Mosquito population-dynamics parameters
    MosqLifeCycleParams lcParams;
    
    
    /** @brief Inputs which are constant after simulation start
     * 
     * Probabilities have no units; others have units specified.
     *
     * All parameters are calculated during initialisation and in theory don't
     * need checkpointing. */
    //@{
    /** Death rate of mosquitoes while host-seeking (μ_vA).
     *
     * TODO: the model could be extended to allow this and mosqSeekingDuration to vary over the year.
     *
     * Unit: animals/day. */
    double mosqSeekingDeathRate;

    /** Probability of a mosquito successfully laying eggs given that it has
     * rested (P_E).
     *
     * Currently assumed constant, although NC's non-autonomous model provides
     * an alternative. */
    double probMosqSurvivalOvipositing;

    struct NHHParams {
        // α_i
        // rate: humans encountered per day
        double entoAvailability;
        // α_i * P_B_i * P_C_i * P_D_i
        // units as for entoAvailability
        double probCompleteCycle;
    };
    /** Non-human host data. Doesn't need checkpointing. */
    vector<NHHParams> nonHumans;

    /// If less than this many mosquitoes remain infected, transmission is interrupted.
    double minInfectedThreshold;
    //@}

    
    ///@brief Descriptions of transmission, used primarily during warmup
    //@{
    /// Angle (in radians) to rotate series generated by FSCoeffic by, for EIR.
    double EIRRotateAngle;

    /// Rotation angle (in radians) for emergence rate. Both offset for EIR given in XML file and
    /// offset needed to fit target EIR (delayed from emergence rate). Checkpoint.
    double FSRotateAngle;

    /** Fourier coefficients for EIR / forcedS_v series, input from XML file.
     *
     * Initially used to calculate initialisation EIR, then scaled to calc. S_v.
     *
     * When calcFourierEIR is used to produce an EIR from this over 365
     * (365) elements, the resulting EIR has units of
     * infectious bites per adult per day.
     *
     * fcEir must have odd length and is ordered: [a0, a1, b1, ..., an, bn].
     * FSCoeffic[0] needs checkpointing, the rest doesn't. */
    vector<double> FSCoeffic;

    /** S_v used to force an EIR during vector init.
     * 
     * Has annual periodicity: length is 365. First value (index 0) corresponds
     * to first day of year (1st Jan or something else if rebased). In 5-day
     * time-step model values at indecies 0 through 4 are used to calculate the
     * state at time-step 1.
     *
     * Should be checkpointed. */
    vector<double> forcedS_v;

    /** Summary of S_v over the last five years, used by vectorInitIterate to
     * calculate scaling factor.
     *
     * Length is 365 * 5. Checkpoint.
     *
     * Units: inoculations. */
    vector<double> quinquennialS_v;

    /** Conversion factor from forcedS_v to mosqEmergeRate.
     *
     * Also has another temporary use between initialise and setupNv0 calls:
     * "initOvFromSv" or  (ρ_O / ρ_S).
     *
     * Should be checkpointed. */
    double initNv0FromSv;       ///< ditto

    /** Conversion factor from forcedS_v to (initial values of) N_v (1 / ρ_S).
     * Should be checkpointed. */
    double initNvFromSv;
    //@}
    
#if 0
    NOTE: we don't want this anymore?
    /** Emergence rate of new mosquitoes, for every day of the year (N_v0).
     * larvalResources is fitted to this.
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
#endif

    /** @brief Parameter arrays N_v_length long.
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

    
    ///@brief Other variables storing state of model
    //@{
    MosquitoLifeCycle lcModel;
    
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
    vector<double> fArray;
    vector<double> ftauArray;
    
    /** Per time-step partial calculation of EIR.
    *
    * See comment in advancePeriod() for details of how the EIR is calculated.
    *
    * Doesn't need to be checkpointed (is recalculated each step). */
    double partialEIR;
    //@}
    

    /** @brief Intervention parameters
     *
     * Would need to be checkpointed for main simulation; not used during
     * initialisation period (so could be reinitialised). */
    //@{
    /** Timestep at which larviciding effects dissappear. */
    TimeStep larvicidingEndStep;
    /** One-minus larviciding effectiveness. I.e. emergence rate is multiplied by
     * this parameter. */
    double larvicidingIneffectiveness;
    //@}
    
    ///@brief Storage for summary data
    /** Variables tracking data to be reported. */
    double timestep_N_v0;
    
    friend class VectorEmergenceSuite;
    friend class SpeciesModelSuite;
};

}
}
}
#endif
