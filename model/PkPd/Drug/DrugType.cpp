/*

  This file is part of OpenMalaria.
 
  Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "PkPd/Drug/DrugType.h"
#include "util/errors.hpp"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

namespace OM { namespace PkPd {
    
/*
 * Static variables and functions
 */

map<string,DrugType*> DrugType::available; 



void DrugType::init () {
}


void DrugType::addDrug(DrugType* drug) {
  string abbrev = drug->abbreviation;
  // Check drug doesn't already exist
    if (available.find (abbrev) != available.end())
    throw invalid_argument (string ("Drug already in registry: ").append(abbrev));
  
  available.insert (pair<string,DrugType*>(abbrev, drug));
}

const DrugType* DrugType::getDrug(string _abbreviation) {
  map<string,DrugType*>::const_iterator i = available.find (_abbreviation);
  if (i == available.end())
    throw util::xml_scenario_error (string ("prescribed non-existant drug ").append(_abbreviation));
  
  return i->second;
}


// -----  Non-static DrugType functions  -----

DrugType::DrugType (string _abbreviation) {
  abbreviation = _abbreviation;
}
DrugType::~DrugType () {}

} }