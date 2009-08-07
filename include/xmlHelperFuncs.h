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

#ifndef Hmod_XmlHelperFuncs
#define Hmod_XmlHelperFuncs

#include "global.h"
#include "scenario.hxx"
#include <vector>

/** Reads a generic list of doubles from XML.
 *
 * @param list XML element to read from.
 * @param length Expected length of list. Will throw if not correct. */
vector<double> readDoubleList (const scnXml::DoubleList& list, size_t length) {
  const scnXml::DoubleList::ItemSequence seq = list.getItem();
  if (seq.size() != length)
    throw xml_scenario_error ("readDoubleList: XML list has wrong length");
  vector<double> ret (length);
  for (size_t i = 0; i < length; ++i)
    ret[i] = seq[i];
  return ret;
}

#endif
