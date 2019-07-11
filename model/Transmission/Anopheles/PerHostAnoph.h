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

#ifndef Hmod_Anopheles_PerHost
#define Hmod_Anopheles_PerHost

#include "util/sampler.h"
#include "schema/entomology.h"

namespace OM {
namespace Transmission {
namespace Anopheles {

/** Stores vector model data applicable between a category of host and a
 * mosquito species: intervention descriptions and model parameters.
 *
 * Parameters are read from XML, and the availability rate is adjusted. */
class PerHostAnophParams {
public:
    static inline void initReserve (size_t numSpecies) {
        params.reserve (numSpecies);
    }
    static inline void init (const scnXml::Mosq& mosq) {
        params.push_back(PerHostAnophParams{ mosq });
    }
    
    /// Get the number of vector species
    static inline size_t numSpecies() {
        return params.size();
    }
    
    /// Get parameters for the given vector species
    static inline const PerHostAnophParams& get(size_t species) {
        return params[species];
    }
    
    /** entoAvailability is calculated externally, then set after other
     * parameters have been initialised.
     * 
     * This function doesn't need to exist, but helps make this fact obvious.
     * 
     * It should be called exactly once.
      */
    inline static void scaleEntoAvailability(size_t species, double entoAvailability){
        params[species].entoAvailability.scaleMean( entoAvailability );
    }

    /** @brief Probabilities of finding a host and surviving a feeding cycle
     * 
     * These parameters describe the mean and heterogeneity of α_i, P_B_i,
     * P_C_i and P_D_i across the human population. */
    //@{
    /** Availability rate (α_i) */
    util::LognormalSampler entoAvailability;

    /** Probability of mosquito successfully biting host (P_B_i) */
    util::BetaSampler probMosqBiting;

    /** Probability of mosquito escaping human and finding a resting site without
     * dying, after biting the human (P_C_i). */
    util::BetaSampler probMosqFindRestSite;

    /** Probability of mosquito successfully resting after finding a resting site
     * (P_D_i). */
    util::BetaSampler probMosqSurvivalResting;
    //@}
    
private:
    PerHostAnophParams (const scnXml::Mosq& mosq);

    static vector<PerHostAnophParams> params;
};

/** Data needed for each human which is per-mosquito species. */
class PerHostAnoph
{
public:
    /** In lieu of a constructor initialises elements, using the passed base to
     * get baseline parameters. */
    void initialise (size_t species, double availabilityFactor);
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        entoAvailability & stream;
        probMosqBiting & stream;
        probMosqRest & stream;
    }
    
    /** Return the availability rate (α_i) of this human to mosquitoes. */
    inline double getEntoAvailability() const {
        return entoAvailability;
    }
    /** Return the probability of mosquito successfully biting host (P_B_i) */
    inline double getProbMosqBiting() const {
        return probMosqBiting;
    }
    /** Return the probability of mosquito escaping human and finding a resting
     * site, then resting without dying, after biting the human (P_C_i * P_D_i). */
    inline double getProbMosqRest() const {
        return probMosqRest;
    }
    
private:
    ///@brief Rate/probabilities before interventions. See functions.
    //@{
    /** Availability rate of human to mosquitoes, including hetergeneity factor
     * and base rate, but excluding age and intervention factors. */
    double entoAvailability;
    
    /** Probability of mosquito successfully biting host (P_B_i) in the absense of
     * interventions. */
    double probMosqBiting;
    
    /** Probability of mosquito escaping human and finding a resting site, then
     * resting without dying, after biting the human (P_C_i * P_D_i) in the
     * absense of interventions. */
    double probMosqRest;
    //@}
};


}
}
}
#endif
