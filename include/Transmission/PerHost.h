/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#ifndef Hmod_PerHost
#define Hmod_PerHost

#include "Transmission/Anopheles/PerHost.h"
#include "Transmission/ITN.h"
#include "Transmission/IRS.h"
#include "util/AgeGroupInterpolation.h"
#include "util/DecayFunction.h"
#include <boost/shared_ptr.hpp>

namespace OM {
namespace Transmission {
    
class TransmissionModel;
using util::AgeGroupInterpolation;
using util::DecayFunction;
using util::DecayFuncHet;
using boost::shared_ptr;

/** Contains TransmissionModel parameters which need to be stored per host.
 *
 * Currently many members are public and directly accessed. */
class PerHost
{
public:
    /// @brief Static member functions
    //@{
    /** Static initialisation. */
    static void init ();
    /** Static cleanup. */
    static void cleanup ();
    
    static void setVADescription (const scnXml::VectorDeterrent& elt);
    //@}
    
    ///@brief Initialisation / checkpionting
    //@{
    PerHost (const Transmission::TransmissionModel& tm);
    void initialise (TransmissionModel& tm, double availabilityFactor);
    //@}
    
    /// Call once per timestep. Updates net holes.
    inline void update(const ITNParams& params) {
        net.update(params);
    }
    
    ///@brief Intervention controls
    //@{
    inline void removeFromTransmission (bool s){
        outsideTransmission = s;
    }
  
    /// Give individual a new ITN as of time timeStep.
    void setupITN (const TransmissionModel& tm);
    /// Give individual a new IRS as of time timeStep.
    void setupIRS (const TransmissionModel& tm);
    /// Give individual a new VA intervention as of time timeStep.
    void setupVA ();
    
    /// Is individual protected by a VA?
    inline bool hasVAProtection(TimeStep maxInterventionAge)const{
        return timestepVA + maxInterventionAge > TimeStep::simulation;
    }
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
    /// Get a reference to the net
    inline const ITN& getITN() const{
        return net;
    }
    /// Get a reference to the IRS
    inline const IRS& getIRS() const{
        return irs;
    }
    
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
        timestepVA & stream;
        net & stream;
        irs & stream;
    }
    //@}
    
private:
    vector<Anopheles::PerHost> species;
    
    // Determines whether human is outside transmission
    bool outsideTransmission;
    
    // Heterogeneity factor in availability; this is already multiplied into the
    // entoAvailability param stored in HostMosquitoInteraction.
    double _relativeAvailabilityHet;
    
    // (TimeStep::simulation - timestepXXX) is the age of the intervention in
    // the new time-step (that being updated).
    // timestepXXX = TIMESTEP_NEVER means intervention has not been deployed.
    TimeStep timestepVA;
    
    DecayFuncHet hetSampleVA;

    ITN net;
    IRS irs;
    
    static AgeGroupInterpolation* relAvailAge;
    
    // descriptions of decay of interventions
    // set if specific intervention is used
    static shared_ptr<DecayFunction> VADecay;
};

}
}
#endif
