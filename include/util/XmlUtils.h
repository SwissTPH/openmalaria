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

#ifndef Hmod_XmlUtils
#define Hmod_XmlUtils

#include "Global.h"
#include "scenario.hxx"

namespace OM { namespace util { namespace XmlUtils {
  ///@brief Convertions between scnXml::DoubleList and std::vector<double>
  //@{
  
  /** Reads a generic list of doubles from XML.
   *
   * @param list XML element to read from.
   * @param length Expected length of list. Will throw if not correct. */
  vector<double> DoubleList2std (const scnXml::DoubleList& list, size_t length);
  //@}
  
  ///@brief Conversions between XML elements and std::map
  //@{
    /** Read an XML element of by-age-group values into a map.
     *
     * @param outp Output map. Keys are lower age bound, lowest is guaranteed
     * to be 0.
     * @param ageGroups XML element to read from.
     * @param eltName Name of XML element; used for error reporting.
     * @param addUbound If true, a final element with age-bound +inf and
     * value equal to that of the largest age group is added to outp.
     */
    void lboundGroups2map (
	map<double,double>& outp,
	const scnXml::AgeGroupValues::GroupSequence& ageGroups,
	const char* eltName,
	bool addUbound = false
    );
  //@}

} } }
#endif
