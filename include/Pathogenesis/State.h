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

#ifndef MO_PATHOGENISIS_STATE
#define MO_PATHOGENISIS_STATE

namespace OM { namespace Pathogenesis {
  /** Types of sickness; used by case management.
   *
   * Most values are flags which can be combined in any form. A few
   * combinations set follow. */
  enum State {
    /* Values here are written in hexadecimal: http://en.wikipedia.org/wiki/Hexadecimal
     * Many are designed to be "flags", so the value corresponds to a single bit:
     * http://en.wikipedia.org/wiki/Flag_byte
     * Max: 0x4000
     * (note & | ^ are C++'s binary AND, OR and XOR operators). */
    NONE		= 0,		///< Not sick
    
    // Flags for current state/worst state to report:
    SICK		= 0x1,		///< Sick (may or may not be from malaria)
    MALARIA		= 0x2,		///< Malaria sickness
    SEVERE		= 0x8,		///< Severe malaria case
    COINFECTION		= 0x4,		///< Malaria with a coinfection
    /// Used by ClinicalEventScheduler to indicate a second bout of malarial sickness within the same episode (roughly)
    SECOND_CASE		= 0x10,
    COMPLICATED		= 0x20,		///< Flag used to indicate SEVERE and/or COINFECTION
    
    MORBIDITY_MASK	= 0x3F,		///< Mask coving all above states
    
    // Flag used by pathogenesis model to tell the clinical model that individual will die; not used for reporting:
    INDIRECT_MORTALITY	= 0x800,	///< Death caused by indirect effects of malaria
    
    // Flags for outcome reporting:
    EVENT_IN_HOSPITAL	= 0x400,	///< Indicates recovery/sequelae/death event occurred in hospital − only set on one of these events (ImmediateOutcomes only)
    DIRECT_DEATH	= 0x1000,	///< Used for reporting death (from COMPLICATED sickness)
    SEQUELAE		= 0x2000,	///< Reporting recovered with sequelae (from COMPLICATED sickness)
    RECOVERY		= 0x4000,	///< Report that individual fully recovered
    
    STATE_MALARIA	= SICK | MALARIA,	///< Combination: SICK, MALARIA
    STATE_SEVERE	= STATE_MALARIA | COMPLICATED | SEVERE,	///< Combination: SICK, MALARIA, COMPLICATED, SEVERE
    STATE_COINFECTION	= STATE_MALARIA | COMPLICATED | COINFECTION,	///< Combination: SICK, MALARIA, COMPLICATED, COINFECTION
  };

} }
#endif
