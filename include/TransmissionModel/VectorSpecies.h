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

#ifndef Hmod_VectorTransmissionSpecies
#define Hmod_VectorTransmissionSpecies

namespace scnXml {
  class Anopheles;
}

/** Per-species data for vector control.
 *
 * A list of this class is used by the VectorTransmission class to hold
 * (potentially) species-specific per-population data.
 * 
 * Instead of storing static variables in this class, store them in
 * the VectorTransmission class. */
class VectorTransmissionSpecies
{
  friend class VectorTransmission;
private:
  /** Called to initialise variables instead of a constructor. */
  void setAnophelesData (scnXml::Anopheles anoph);
  
  /** @brief Parameters describing the mosquitos' success in NC's model.
   *
   * Apart from the ovipositing probability, these may be varied per-human to
   * account for interventions and innate resistances. */
  //@{
  //FIXME: not set:
  double entoAvailability;
  
  /** Probability of mosquito successfully biting host (P_B_i) */
  double probMosqBiting;
  
  /** Probability of mosquito escaping human and finding a resting site without
   * dying, after biting the human (P_C_i). */
  double probMosqFindRestSite;
  
  /** Probability of mosquito successfully resting after finding a resting site
   (P_D_i). */
  double probMosqSurvivalResting;
  
  /** Probability of a mosquito successfully laying eggs given that it has
  * rested (P_E).
  * 
  * Currently assumed constant, although NC's non-autonomous model provides
  * an alternative. */
  double probMosqSurvivalOvipositing;
  //@}
};

#endif
