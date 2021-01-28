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

#ifndef Hmod_AnophelesModel
#define Hmod_AnophelesModel

#include "Global.h"
#include "Transmission/Anopheles/PerHostAnoph.h"
#include "Transmission/PerHost.h"
#include "util/SimpleDecayingValue.h"
#include "util/vectors.h"
#include "util/vecDay.h"

#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

#include "rotate.h"

#include <vector>
#include <limits>

namespace OM {
namespace Transmission {
namespace Anopheles {
    using std::numeric_limits;
    using util::vector2D;
    using namespace OM::util;

enum VecStat { PA, PAmu, PA1, PAh, PDF, PDIF, NV, OV, SV };

/** Per-species part for vector transmission model.
 *
 * Data in this class is specific to a species of anopheles mosquito, where
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
 * a heterogeneous host population" (3rd Oct. 2007).
 * 
 * Some terminology extensions (not in paper):
 * 
 * $\nu_A$ / ν_A / leaveRate = sum_i=1^n α_i N_i + μ_vA
 * $\alpha_d$ / α_d / availDivisor = P_Ai / (α_i N_i) = (1 - P_A) / ν_A
 * $\sigma_{df}$ / σ_f / sigma_f = sum_i α_i N_i P_Bi
 * $\sigma_{df}$ / σ_df / sigma_df = sum_i α_i N_i P_Bi P_Ci P_Di
 * $\sigma_{dif}$ / σ_dif / sigma_dif = sum_i α_i N_i P_Bi P_Ci P_Di K_vi
 * 
 * Thus:
 * 
 * P_df = σ_df α_d P_E
 * P_dif = σ_dif α_d P_E
 */

struct MosquitoParams
{
    // Proportion of host-seeking parous mosquitoes (those who have laid eggs) which laid eggs that day:
    double laidEggsSameDayProportion;
    // Probability that the mosquito survives the feeding cycle.
    // Note: Pf = M, the parous rate (prop mosqs which have laid eggs):
    double survivalFeedingCycleProbability;
    double humanBloodIndex;    // χ (chi)
    // Cycle probabilities, when biting a human:
    // B: Host encountered
    double probBiting;
    // C: Fed
    double probFindRestSite;
    // D: Resting
    double probResting;
    // E: Laying eggs (ovipositing)
    double probOvipositing;

    /** Duration of host-seeking per day; the maximum fraction of a day that a
     * mosquito would spend seeking (θ_d). */
    double seekingDuration;

    /** Probability of a mosquito successfully laying eggs given that it has
     * rested (P_E).
     *
     * Currently assumed constant, although NC's non-autonomous model provides
     * an alternative. */
    double probMosqSurvivalOvipositing;

    /** If less than this many mosquitoes remain infected, transmission is interrupted. */
    double minInfectedThreshold;

    /** @brief Probabilities and rates associated with life-cycle model
     * 
     * These are calculated during initialisation and thereafter constant.
     * 
     * Probabilities have no units; others have units specified.
     *
     * All parameters are calculated during initialisation and in theory don't
     * need checkpointing. 
     * Death rate of mosquitoes while host-seeking (μ_vA).
     *
     * Unit: animals/day. */
    double seekingDeathRate;

    /** Duration parameters for mosquito/parasite life-cycle
     * 
     * Currently these are all constant. In theory they could be made to vary
     * seasonally, based on a fixed periodic cycle, though some code and
     * possibly model changes would be needed to accomodate this.
     * 
     * All have units of days.
     *
     * Set in initialise function from XML data; no need to checkpoint.
     * Duration of feeding cycle (equals duration of resting period) for mosquito (τ).
     * Units: days. */
    SimTime restDuration;

    /** Duration of the extrinsic incubation period (sporozoite development time)
    * (θ_s).
    * Units: Days.
    *
    * Doesn't need checkpointing. */
    SimTime EIPDuration;

    string name;
};

struct NhhParams {
    double mosqRelativeEntoAvailability;
    double mosqProbBiting;
    double mosqProbFindRestSite;
    double mosqProbResting;
    double hostFecundityFactor;
    string name;
};

struct NhhParamsInterv {
    double mosqRelativeAvailabilityHuman;
    double mosqProbBiting;
    double mosqProbFindingRestSite;
    double mosqProbResting;
    double hostFecundityFactor;
};

struct Nhh {
    double avail_i;
    double P_B_I;
    double P_C_I;
    double P_D_I;
    double rel_fecundity;
    SimTime expiry;
};

struct TrapParams {
    TrapParams(): relAvail(numeric_limits<double>::signaling_NaN()) {}
    TrapParams(TrapParams&& o): relAvail(o.relAvail), availDecay(move(o.availDecay)) {}
    
    double relAvail; // Initial availability of a trap relative to an adult
    unique_ptr<util::DecayFunction> availDecay; // Decay of availability
};

struct TrapData {
    size_t instance; // index in trapParams
    double initialAvail; // initial availability (avail per trap * num traps)
    DecayFuncHet availHet; // parameter for decay of availability
    SimTime deployTime; // deploy time (for decay function)
    SimTime expiry; // date at which this intervention should be deleted
};

class AnophelesModel
{
public:
    ///@brief Initialisation and destruction
    //@{
    AnophelesModel () :
            nhh_avail(numeric_limits<double>::signaling_NaN()),
            nhh_sigma_df(numeric_limits<double>::signaling_NaN()),
            nhh_sigma_dff(numeric_limits<double>::signaling_NaN()),
            // EmergenceModel
            EIRRotateAngle(numeric_limits<double>::quiet_NaN()),
            initNvFromSv(numeric_limits<double>::quiet_NaN()),
            initOvFromSv(numeric_limits<double>::quiet_NaN()),
            initNv0FromSv(numeric_limits<double>::quiet_NaN()),
            interventionSurvival(1.0),
            // MosqTransmission
            N_v_length(SimTime::zero()),
            timeStep_N_v0(0.0)
    {
        forcedS_v.resize (SimTime::oneYear());
        quinquennialS_v.assign (SimTime::fromYearsI(5), 0.0);
        mosqEmergeRate.assign (SimTime::oneYear(), 0.0);
    }
    
    AnophelesModel(const AnophelesModel&) = delete;            //disable copy-constructor
    AnophelesModel& operator=(const AnophelesModel&) = delete; //disable copy-assignment
    AnophelesModel(AnophelesModel&&) = default;

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
    void initialise ( size_t species, MosquitoParams mosqParams);
    
    ///@brief Initialisation helper functions
    //@{
    /** Calculate availability rate of hosts (α_i) and death rate while seeking
     * (µ_vA)
     *
     * Documentation: "Parameter Values for Transmission model"
     * (Chitnis, Smith and Schapira, 4.3.2010)
     * 
     * @param species Species index
     * @param anoph Data from XML
     * @param nonHumanHostPopulations Size of each non-human population
     * @param populationSize Size of the human population (assumed constant)
     */
    void initAvailability(size_t species, const vector<NhhParams> &nhhs, int populationSize);

    void initEIR(vector<double>& initialisationEIR, vector<double> FSCoefficInit, double EIRRotateAngleInit, double targetEIR, double propInfectious, double propInfected);

    /** Scale the internal EIR representation by factor; used as part of
     * initialisation. */
    inline void scaleEIR( double factor ){
        FSCoeffic[0] += log( factor );
    }

    virtual void scale(double factor)
    {
        vectors::scale (N_v, factor);
        vectors::scale (O_v, factor);
        vectors::scale (S_v, factor);
    }
    
    /** Initialisation which must wait until a human population is available.
     * This is only called when a checkpoint is not loaded.
     *
     * @param nHumans Human population size
     * @param meanPopAvail The mean availability of age-based relative
     * availability of humans to mosquitoes across populations.
     * @param sum_avail sum_i α_i * N_i for human hosts i
     * @param sigma_f sum_i α_i * N_i * P_Bi for human hosts i
     * @param sigma_df sum_i α_i * N_i * P_Bi * P_Ci * P_Di for human hosts i
     * @param sigma_dff sum_i α_i * N_i * P_Bi * P_Ci * P_Di * rel_mosq_fecundity for human hosts i
     *
     * Can only usefully run its calculations when not checkpointing, due to
     * population not being the same when loaded from a checkpoint. */
    virtual void init2 (int nHumans, double meanPopAvail,
        double sum_avail, double sigma_f, double sigma_df, double sigma_dff);
    
    /** Set up the non-host-specific interventions. */
    void initVectorInterv( const scnXml::VectorSpeciesIntervention& elt, size_t instance );
    
    /** Set up trap parameters. */
    void initVectorTrap( const scnXml::Description1& desc, size_t instance );
    
    // /** Set up the non-human hosts interventions. */
    void initNonHumanHostsInterv( const scnXml::NonHumanHostsSpeciesIntervention& elt, const scnXml::DecayFunction& decay, size_t instance, string name );
    
    void initAddNonHumanHostsInterv( const scnXml::NonHumanHostsVectorSpecies& elt, string name );
    
    /** Work out whether another interation is needed for initialisation and if
     * so, make necessary changes.
     *
     * @returns true if another iteration is needed. */
    virtual bool initIterate();

    ///@brief Functions called as part of usual per-time-step operations
    //@{
    /** Called per time-step. Does most of calculation of EIR.
     *
     * @param sum_avail sum_i α_i * N_i for human hosts i
     * @param sigma_df sum_i α_i * N_i * P_Bi * P_Ci * P_Di for human hosts i
     * @param sigma_dif sum_i α_i * N_i * P_Bi * P_Ci * P_Di * Kvi for human hosts i; this vector is modified!
     * @param sigma_dff sum_i α_i * N_i * P_Bi * P_Ci * P_Di * rel_mosq_fecundity for human hosts i
     * @param isDynamic True to use full model; false to drive model from current contents of S_v.
     */
    void advancePeriod (double sum_avail, double sigma_df, vector<double>& sigma_dif, double sigma_dff, bool isDynamic);

    /// Intermediatary from vector model equations used to calculate EIR
    inline const vector<double>& getPartialEIR() const{ return partialEIR; }
    //@}


    ///@brief Functions called to deploy interventions
    //@{
    void deployVectorPopInterv (LocalRng& rng, size_t instance);
    /// Deploy some traps
    /// 
    /// @param instance Index of this type of trap
    /// @param number The number of traps to deploy
    /// @param lifespan Time until these traps are removed/replaced/useless
    void deployVectorTrap(LocalRng& rng, size_t species, size_t instance, double popSize, SimTime lifespan);

    void deployNonHumanHostsInterv(LocalRng& rng, size_t species, size_t instance, string name);

    void deployAddNonHumanHosts(LocalRng& rng, size_t species, string name, double popSize, SimTime lifespan);
    
    /** (Re) allocate and initialise some state variables. Must be called
     * before model is run. */
    void initState ( double tsP_A, double tsP_Amu, double tsP_A1, double tsP_Ah,
                     double tsP_df, double tsP_dff,
                     double initNvFromSv, double initOvFromSv,
                     const vecDay<double>& forcedS_v);
    
    /// Helper function for initialisation.
    void initIterateScale ( double factor );
    
    virtual double getEmergenceRate(const SimTime &d0, const vecDay<double>& mosqEmergeRate, double nOvipositing)
    {   
        // Get emergence at start of step:
        SimTime dYear1 = mod_nn(d0, SimTime::oneYear());
        // Simple model: fixed emergence scaled by larviciding
        return mosqEmergeRate[dYear1];
    }

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
    void update( SimTime d0, double tsP_A, double tsP_Amu, double tsP_A1, double tsP_Ah, double tsP_df,
                   const vector<double> tsP_dif, double tsP_dff,
                   bool isDynamic,
                   vector<double>& partialEIR, double EIR_factor);
    
    ///@brief Interventions and reporting
    //@{
    void uninfectVectors();
    //@}
    
    inline SimTime getEIPDuration() const {
        return mosq.EIPDuration;
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
    
    virtual double getResAvailability() const {
        return numeric_limits<double>::quiet_NaN();
    }
    virtual double getResRequirements() const {
        return numeric_limits<double>::quiet_NaN();
    }
    
    /// Write some per-species summary information.
    void summarize( size_t species )const;
    //@}


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
    //@}
    
    virtual void checkpoint (istream& stream){ (*this) & stream; }
    virtual void checkpoint (ostream& stream){ (*this) & stream; }

    /// Checkpointing
    //Note: below comments about what does and doesn't need checkpointing are ignored here.
    template<class S>
    void operator& (S& stream) {
        mosq.seekingDeathRate & stream;
        mosq.seekingDuration & stream;
        mosq.probMosqSurvivalOvipositing & stream;
        // transmission & stream;
        seekingDeathRateIntervs & stream;
        probDeathOvipositingIntervs & stream;
        partialEIR & stream;
        EIRRotateAngle & stream;
        FSCoeffic & stream;
        forcedS_v & stream;
        initNvFromSv & stream;
        initOvFromSv & stream;
        emergenceReduction & stream;
        interventionSurvival & stream;
        // EmergenceModel
        mosqEmergeRate & stream;
        quinquennialS_v & stream;
        initNv0FromSv & stream;
        // (*emergence) & stream;
        // MosqTransmission
        mosq.restDuration & stream;
        mosq.EIPDuration & stream;
        N_v_length & stream;
        P_A & stream;
        P_df & stream;
        P_dif & stream;
        P_dff & stream;
        N_v & stream;
        O_v & stream;
        S_v & stream;
        P_Amu & stream;
        P_A1 & stream;
        P_Ah & stream;
        //TODO: do we actually need to checkpoint these next three?
        fArray & stream;
        ftauArray & stream;
        uninfected_v & stream;
        timeStep_N_v0 & stream;
    }

    MosquitoParams mosq;

    // sum_i N_i * α_i for i in NHH types: total availability of non-human hosts
    double nhh_avail;
    // sum_i N_i * α_i * P_B_i * P_C_i * P_D_i for i in NHH types: chance feeding
    // and surviving a complete cycle across all non-human hosts
    double nhh_sigma_df;
    // as above, but modified by fertility factors
    double nhh_sigma_dff;
    //@}
    
    
    // Parameters for trap interventions. Doesn't need checkpointing.
    vector<TrapParams> trapParams;
    
    // Added NonHumanHosts
    map<string,Nhh> nhhInstances;

    map<string,NhhParamsInterv> nhhDefinitionsInterv;

    // -----  model state (and some encapsulated parameters)  -----
    /** @brief Intervention parameters */
    //@{
    /** Interventions affecting death rate while seeking (parameters + state)
     * 
     * Value is increase in death rate (so multiply rate by 1 + this). */
    vector<util::SimpleDecayingValue> seekingDeathRateIntervs;
    /** Interventions affecting probability of dying while ovipositing (pre laying of eggs) (parameters + state)
     * 
     * Value is probability of dying due to this intervention (so multiply survival by 1 - this). */
    vector<util::SimpleDecayingValue> probDeathOvipositingIntervs;

    /** Baited trap interventions.
     * 
     * This is the "full" availability: availability per trap times the number
     * of traps. Each deployment has its own availiability along with expiry
     * date; total availability is the sum. */
    list<TrapData> baitedTraps;
    //@}

    map<string,vector<util::SimpleDecayingValue>> reduceNhhAvailability, reduceP_B_I, reduceP_C_I, reduceP_D_I, reduceFecundity;
    
    /** Per time-step partial calculation of EIR, per genotype.
    *
    * See comment in advancePeriod() for details of how the EIR is calculated.
    *
    * Doesn't need to be checkpointed (is recalculated each step). */
    vector<double> partialEIR;


    // Emergence Model
    // -----  parameters (constant after initialisation)  -----
    
    ///@brief Descriptions of transmission, used primarily during warmup
    //@{
    /// Angle (in radians) to rotate series generated by FSCoeffic by, for EIR.
    double EIRRotateAngle;

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
    vecDay<double> forcedS_v;
    
    /** Conversion factor from forcedS_v to (initial values of) N_v (1 / ρ_S).
     * Should be checkpointed. */
    double initNvFromSv;
    
    /** Conversion factor from forcedS_v to (initial values of) O_v (ρ_O / ρ_S).
     * Should be checkpointed. */
    double initOvFromSv;
    //@}

    ///@brief Descriptions of transmission, used primarily during warmup
    //@{
    /** Summary of S_v over the last five years, used by vectorInitIterate to
     * calculate scaling factor.
     *
     * Length is 365 * 5. Checkpoint.
     *
     * Units: inoculations. */
    vecDay<double> quinquennialS_v;
    
    /** Conversion factor from forcedS_v to mosqEmergeRate.
     *
     * Should be checkpointed. */
    double initNv0FromSv;       ///< ditto
    //@}
    
    /** @brief Intervention parameters
     *
     * Checkpointed. */
    //@{
    /// Description of intervention killing effects on emerging pupae
    vector<util::SimpleDecayingValue> emergenceReduction;
    /// Cache parameter updated by update()
    double interventionSurvival;   // survival with regards to intervention effects
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

    /** Probability of a mosquito dying */
    vecDay<double> P_Amu;

    /** Probability of a mosquito finding a host */
    vecDay<double> P_A1;

    /** Probability of a mosquito finding a non-human host type*/
    vecDay<double> P_Ah;
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
};

}
}
}
#endif
