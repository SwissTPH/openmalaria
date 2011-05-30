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

#ifndef Hmod_HostCategoryAnopheles
#define Hmod_HostCategoryAnopheles

#include "Global.h"
#include "Transmission/Vector/ITN.h"

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
    probMosqFindRestSite(0.0), probMosqSurvivalResting(0.0)
  {}
  
  ~HostCategoryAnopheles() {}
  
  /** set the (human or non human) ento Availability.
   *  This is only a helper method, since the ento availability
   *  is calculated in VectorAnopheles.
   */
  void setEntoAvailability(double entoAvailability);
  
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
};


} }
#endif
