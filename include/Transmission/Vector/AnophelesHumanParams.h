/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef HOSTCATEGORYANOPHELESHUMANS_H_
#define HOSTCATEGORYANOPHELESHUMANS_H_

#include "Transmission/Vector/ITN.h"
#include "util/sampler.h"

namespace OM {
namespace Transmission {

/** Stores vector model data applicable between a category of host and a
 * mosquito species: intervention descriptions and model parameters.
 *
 * Parameters are read from XML, and the availability rate is adjusted. */
class AnophelesHumanParams {
public:
    AnophelesHumanParams(const ITNParams* baseITNParams) :
            net( baseITNParams ),
            IRSDeterrency(numeric_limits< double >::signaling_NaN()),
            IRSKillingEffect(numeric_limits< double >::signaling_NaN()),
            VADeterrency(numeric_limits< double >::signaling_NaN())
    {}
    
    /** Set parameters from an XML element. */
    void operator= (const scnXml::Mosq& mosq);

    /** entoAvailability is calculated externally, then set after other
     * parameters have been initialised.
     * 
     * This function doesn't need to exist, but helps make this fact obvious.
      */
    inline void setEntoAvailability(double entoAvailability){
        this->entoAvailability.setMean( entoAvailability );
    }
    
    /** Set up vector-model intervention parameters. */
    void setITNDescription (const ITNParams& params, const scnXml::ITNDescription::AnophelesParamsType& elt, double proportionUse);
    /** Set up vector-model intervention parameters. */
    void setIRSDescription (const scnXml::IRSDescription& irsDesc);
    /** Set up vector-model intervention parameters. */
    void setVADescription (const scnXml::BaseInterventionDescription& vaDesc);
    

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
    
    /** @brief Intervention description parameters
     *
     * These describe initial effectiveness. Decay rate/shape is specified
     * elsewhere (by DecayFunction type). */
    //@{
    ITNAnophelesParams net;
    double IRSDeterrency;
    double IRSKillingEffect;
    double VADeterrency;
    //@}
};


}
}
#endif /* HOSTCATEGORYANOPHELESHUMANS_H_ */
