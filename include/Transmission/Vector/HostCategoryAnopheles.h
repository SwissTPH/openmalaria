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

#ifndef Hmod_HostCategoryAnopheles
#define Hmod_HostCategoryAnopheles

#include "Global.h"
#include "inputData.h"

#include <stdexcept>
#include <string>
#include <limits>

namespace OM { namespace Transmission {
    
/** Stores vector model data applicable between a category of host and a
 * mosquito species.
 * 
 * This is a superclass of HostCategoryAnophelesHumans and HostCategoryAnophelesNonHumans
 * which contains parameters found in Human hosts and non human hosts
   *
   * Read from XML by VectorTransmission constructor. No need to checkpoint. */
class HostCategoryAnopheles
{
public:
  /** Initialise entoAvailability and probMosq... to 0. */
  HostCategoryAnopheles () :
    entoAvailability(0.0), probMosqBiting(0.0),
    probMosqFindRestSite(0.0), probMosqSurvivalResting(0.0),
    ITNDeterrency(numeric_limits< double >::signaling_NaN()),
    ITNPreprandialKillingEffect(numeric_limits< double >::signaling_NaN()),
    ITNPostprandialKillingEffect(numeric_limits< double >::signaling_NaN()),
    IRSDeterrency(numeric_limits< double >::signaling_NaN()),
    IRSKillingEffect(numeric_limits< double >::signaling_NaN()),
    VADeterrency(numeric_limits< double >::signaling_NaN())
  {}
  
  ~HostCategoryAnopheles() {}
  
  /** set the (human or non human) ento Availability.
   *  This is only a helper method, since the ento availability
   *  is calculated in VectorAnopheles.
   */
  void setEntoAvailability(double entoAvailability);

    /** Set up vector-model intervention parameters. */
    void setITNDescription (const scnXml::ITNDescription& itnDesc);
    /** Set up vector-model intervention parameters. */
    void setIRSDescription (const scnXml::IRSDescription& irsDesc);
    /** Set up vector-model intervention parameters. */
    void setVADescription (const scnXml::BaseInterventionDescription& vaDesc);
  
  inline double probMosqBitingAndResting() const {
    return probMosqBiting * probMosqFindRestSite * probMosqSurvivalResting;
  }
  
  /// @brief Probabilities of finding a host and surviving a feeding cycle
  //@{
  /** Availability rate (α_i) */
  double entoAvailability;
  
  /** Probability of mosquito successfully biting host (P_B_i) */
  double probMosqBiting;
  
  /** Probability of mosquito escaping human and finding a resting site without
   * dying, after biting the human (P_C_i). */
  double probMosqFindRestSite;
  
  /** Probability of mosquito successfully resting after finding a resting site
   * (P_D_i). */
  double probMosqSurvivalResting;
  //@}

  /** @brief Intervention description parameters
   *
   * These describe initial effectiveness. Decay rate/shape is specified
   * elsewhere (by DecayFunction type). */
  //@{
  /** Effectiveness of net in preventing a mosquito from finding an individual,
   * but not killing the mosquito. (1 - this) multiplies availability. */
  double ITNDeterrency;
  /** (1 - this) is the proportion of mosquitoes killed when trying to feed on
   * an individual. */
  double ITNPreprandialKillingEffect;
  /** (1 - this) is the proportion of mosquitoes killed when trying to escape
   * after feeding on an individual. */
  double ITNPostprandialKillingEffect;
  /** Effectiveness of IRS in preventing a mosquito from finding an individual,
   * but not killing the mosquito. (1 - this) multiplies availability. */
  double IRSDeterrency;
  /** (1 - this) is the proportion of mosquitoes killed when trying to rest. */
  double IRSKillingEffect;
  /** Effectiveness of [intervention] in preventing a mosquito from finding an individual,
   * but not killing the mosquito. (1 - this) multiplies availability. */
  double VADeterrency;
  //@}
};


} }
#endif
