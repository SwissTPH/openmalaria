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

#include "util/WeibullDecayedValue.h"
#include "scenario.hxx"

/** Stores vector model data applicable between a category of host and a
 * mosquito species.
 * 
 * This is either base parameters for a mosquito species to humans, or (final)
 * parameters for a mosquito species to a non-human host type (since non-human
 * hosts are not modelled as agents).
   *
   * Read from XML by VectorTransmission constructor. No need to checkpoint. */
class HostCategoryAnopheles
{
public:
  /** Initialise entoAvailability and probMosq... to 0. */
  HostCategoryAnopheles () :
    entoAvailability(0.0), probMosqBiting(0.0),
    probMosqFindRestSite(0.0), probMosqSurvivalResting(0.0)
  {}
  
  void operator= (const scnXml::Mosq&);
  void operator= (const scnXml::NonHumanHosts&);
  
  /** Set an ITN description for this anopheles species. */
  inline void setITNDescription (const scnXml::Anopheles1& itnDesc) {
    ITNDeterrency = itnDesc.getDeterrency ();
    ITNPreprandialKillingEffect = itnDesc.getPreprandialKillingEffect ();
    ITNPostprandialKillingEffect = itnDesc.getPostprandialKillingEffect ();
  }
  
  /** Set an IRS description for this anopheles species. */
  inline void setIRSDescription (const scnXml::Anopheles2& irsDesc) {
    IRSDeterrency = irsDesc.getDeterrency ();
    IRSKillingEffect = irsDesc.getKillingEffect ();
  }
  
  inline void setVADescription (const scnXml::Anopheles3& vaDesc) {
    VADeterrency = vaDesc.getDeterrency ();
  }
  
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
  
  /** @brief Intervention description parameters */
  //@{
  /** Effectiveness of net in preventing a mosquito from finding an individual,
   * but not killing the mosquito. (1 - this) multiplies availability. */
  WeibullDecayedValue ITNDeterrency;
  /** (1 - this) is the proportion of mosquitoes killed when trying to feed on
   * an individual. */
  WeibullDecayedValue ITNPreprandialKillingEffect;
  /** (1 - this) is the proportion of mosquitoes killed when trying to escape
   * after feeding on an individual. */
  WeibullDecayedValue ITNPostprandialKillingEffect;
  /** Effectiveness of IRS in preventing a mosquito from finding an individual,
   * but not killing the mosquito. (1 - this) multiplies availability. */
  WeibullDecayedValue IRSDeterrency;
  /** (1 - this) is the proportion of mosquitoes killed when trying to rest. */
  WeibullDecayedValue IRSKillingEffect;
  /** Effectiveness of [intervention] in preventing a mosquito from finding an individual,
   * but not killing the mosquito. (1 - this) multiplies availability. */
  WeibullDecayedValue VADeterrency;
  //@}
};
typedef vector<HostCategoryAnopheles> NonHumanHostsType;

#endif
