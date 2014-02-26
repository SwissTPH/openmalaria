/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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
#ifndef Hmod_PerHost
#define Hmod_PerHost

#include "Transmission/Anopheles/PerHost.h"
#include "interventions/Interfaces.hpp"
#include "util/AgeGroupInterpolation.h"
#include "util/DecayFunction.h"
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/shared_ptr.hpp>

namespace OM {
namespace Transmission {
    
class TransmissionModel;
using util::AgeGroupInterpolation;
using util::DecayFunction;
using util::DecayFuncHet;
using boost::shared_ptr;

class HumanVectorInterventionEffect;

/**
 * A base class for interventions affecting human-vector interaction.
 * 
 * The constructor should initialise the data to represent an intervention
 * deployed at this time (TimeStep::simulation).
 * 
 * redeploy() should reset the intervention to a freshly deployed state. If
 * necessary, PerHost::deployEffect can be updated to make it create a new
 * instance instead of calling redeploy.
 */
class PerHostInterventionData {
public:
    /** Deploy an intervention. */
    virtual void redeploy( const HumanVectorInterventionEffect& params ) =0;
    
    /** Per-timestep update. Used by ITNs to update hole decay. */
    virtual void update() =0;
    
    /** Get effect of deterrencies of interventions, as an attractiveness multiplier.
     * 
     * @return a value describing effect on attractiveness. Must not be
     * negative. 0 means mosquitoes are fully deterred, 1 that the intervention
     * has no effect, 2 that the intervention attracts twice as many mosquitoes
     * as would otherwise come. */
    virtual double relativeAttractiveness(size_t speciesIndex) const =0;
    /** Get the killing effect on mosquitoes before they've eaten as a survival
     * multiplier. */
    virtual double preprandialSurvivalFactor(size_t speciesIndex) const =0;
    /** Get the killing effect on mosquitoes after they've eaten as a survival
     * multiplier. */
    virtual double postprandialSurvivalFactor(size_t speciesIndex) const =0;
    
    /// Index of effect describing the intervention
    inline interventions::EffectId id() const { return m_id; }
    
    /// Checkpointing (write only)
    void operator& (ostream& stream) {
        m_id & stream; // must be first; read externally so that the correct
                // makeHumanPart function can be found
        checkpoint( stream );
    }
    
protected:
    /// Set the effect index
    explicit PerHostInterventionData( interventions::EffectId id ) :
            deployTime( TimeStep::simulation ),
            m_id(id) {}
    
    /// Checkpointing: write
    virtual void checkpoint( ostream& stream ) =0;
    
    TimeStep deployTime;        // time of deployment or TimeStep::never
    interventions::EffectId m_id;       // effect index; don't change
};

/** A base class for human vector intervention parameters. */
class HumanVectorInterventionEffect : public interventions::HumanInterventionEffect {
public:
    /** Create a new object to store human-specific details of deployment. */
    virtual PerHostInterventionData* makeHumanPart() const =0;
    virtual PerHostInterventionData* makeHumanPart( istream& stream, interventions::EffectId id ) const =0;
protected:
    explicit HumanVectorInterventionEffect(interventions::EffectId id) : HumanInterventionEffect(id) {}
};

/** Contains TransmissionModel parameters which need to be stored per host.
 *
 * Currently many members are public and directly accessed. */
class PerHost
{
public:
    /// @brief Static member functions
    //@{
    /** Static initialisation. */
    static void init (const scnXml::AgeGroupValues& availabilityToMosquitoes);
    /** Static cleanup. */
    static void cleanup ();
    //@}
    
    ///@brief Initialisation / checkpionting
    //@{
    PerHost (const Transmission::TransmissionModel& tm);
    void initialise (TransmissionModel& tm, double availabilityFactor);
    //@}
    
    /// Call once per timestep. Updates net holes.
    void update();
    
    ///@brief Intervention controls
    //@{
    inline void removeFromTransmission (bool s){
        outsideTransmission = s;
    }
  
    /// Deploy some intervention effect
    void deployEffect( const HumanVectorInterventionEffect& params );
    //@}
    
    /** @brief Availability of host to mosquitoes */
    //@{
    /** Return true if the human has been removed from transmission. */
    inline bool isOutsideTransmission() {
        return outsideTransmission;
    }
    
    /** Calculates the adjustment for body size in exposure to mosquitoes,
     * relative to an average adult.
     * 
     * The bites are assumed proportional to average surface area for hosts of
     * the given age. Linear interpolation is used to calculate this from the
     * input array of surface areas. 
     * 
     * @param ageYears Age of host
     * @return the ratio of bites received by the host to the average for an adult 
     *
     * This is the age factor of availiability; mean output should be
     * mean population availability (that is, 1.0/invMeanPopAvail).
     * 
     * Also has a switch to put individuals entirely outside transmission. */
    inline double relativeAvailabilityAge (double ageYears) const {
        return outsideTransmission ? 0.0 :
            relAvailAge->eval( ageYears );
    }
    
    /** Relative availability of host to mosquitoes excluding age factor.
     *
     * (ONLY for HeterogeneityWorkaroundII, and documentation purposes.)
     * Assume mean is 1.0. */
    inline double relativeAvailabilityHet () const {
        return _relativeAvailabilityHet;
    }
    
    /** Get the availability of this host to mosquitoes relative to an average
     * adult (including heterogeneity and age effects).
     *
     * Used to drive a simulation from an input EIR.
     * Is relativeAvailabilityHet()*relativeAvailabilityAge(ageYears).
     * 
     * Mean output is less than 1.0 (roughly 1.0/invMeanPopAvail).
     */
    inline double relativeAvailabilityHetAge (double ageYears) const {
        return _relativeAvailabilityHet
            * relativeAvailabilityAge (ageYears);
    }
    
    /** Availability of host to mosquitoes (α_i) excluding age factor
     * 
     * (Includes heterogeneity, intervention, and human-to-vector availability
     * rate factors.)
     * 
     * Assume mean is human-to-vector availability rate factor. */
    double entoAvailabilityHetVecItv( const Anopheles::PerHostBase& base, size_t speciesIndex ) const;
    
    /** Availability rate of human to mosquitoes (α_i). Equals 
     * entoAvailabilityHetVecItv()*getRelativeAvailability().
     *
     * To be clear, this includes effects from HetVecItc (het, interv, availability
     * rate) as well as age (avail. relative to an adult). It does not divide
     * by the average availability of the population, which was incorrectly done
     * in the past. */
    inline double entoAvailabilityFull (
        const Anopheles::PerHostBase& base,
        size_t speciesIndex,
        double ageYears
    ) const {
        return entoAvailabilityHetVecItv (base, speciesIndex)
            * relativeAvailabilityAge (ageYears);
    }
    //@}
    
    ///@brief Get killing effects of interventions pre/post biting
    //@{
    /** Probability of a mosquito succesfully biting a host (P_B_i). */
    double probMosqBiting (const Anopheles::PerHostBase& base, size_t speciesIndex) const;
    /** Probability of a mosquito succesfully finding a resting
     * place after biting and then resting (P_C_i * P_D_i). */
    double probMosqResting (const Anopheles::PerHostBase& base, size_t speciesIndex) const;
    /** Set true to remove human from transmission. Must set back to false
     * to restore transmission. */
    //@}
    
    ///@brief Miscellaneous
    //@{
    /** Get the age at which individuals are considered adults (i.e. where
     * availability to mosquitoes reaches its maximum). */
    static inline double adultAge() {
        return relAvailAge->firstGlobalMaximum();
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        species & stream;
        _relativeAvailabilityHet & stream;
        outsideTransmission & stream;
        checkpointIntervs( stream );
    }
    //@}
    
private:
    void checkpointIntervs( ostream& stream );
    void checkpointIntervs( istream& stream );
    
    vector<Anopheles::PerHost> species;
    
    // Determines whether human is outside transmission
    bool outsideTransmission;
    
    // Heterogeneity factor in availability; this is already multiplied into the
    // entoAvailability param stored in HostMosquitoInteraction.
    double _relativeAvailabilityHet;

    typedef boost::ptr_list<PerHostInterventionData> ListActiveEffects;
    ListActiveEffects activeEffects;
    
    static AgeGroupInterpolation* relAvailAge;
};

}
}
#endif
