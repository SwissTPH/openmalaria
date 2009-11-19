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

#include "Drug/Drug.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using namespace std;

/*
 * Static variables and functions
 */

double Drug::minutesPerTimeStep;


void Drug::init () {
  minutesPerTimeStep = Global::interval * 24*60;
}

Drug::Drug(const DrugType* type) :
  typeData (type),
  _concentration (0),
  _nextConcentration (0)
{}
Drug::Drug (const DrugType* type, istream& in) :
  typeData (type)
{
  in >> _concentration;
  // No need to checkpoint; just recalculate _nextConcentration
  _nextConcentration = _concentration * decayFactor (minutesPerTimeStep);
  int num;
  in >> num;
  Global::validateListSize (num);
  for (int i = 0; i < num; ++i)
    doses.push_back (Dose (in));
}

void Drug::write (ostream& out) const {
  out << _concentration << endl;
  assert (doses.size() == 0);	// doses not used yet (remove this eventually)
  out << doses.size() << endl;
  for (deque<Dose>::const_iterator it = doses.begin(); it != doses.end(); ++it)
    it->write (out);
}


void Drug::addDose (double concentration, int delay) {
  if (delay == 0) {
    _concentration += concentration;
    _nextConcentration = _concentration * decayFactor (minutesPerTimeStep);
  } else {
    assert (false);	//NOTE: Not properly dealt with yet (as it is only relevant for ACTs).
    // Only adding doses for this timestep is supported
    assert (delay>0 && delay<minutesPerTimeStep);
    double nextConcentration = concentration*decayFactor (minutesPerTimeStep-delay);
    doses.push_back (Dose (nextConcentration, 0 /*FIXME*/));
  }
}

bool Drug::decay() {
  _concentration = _nextConcentration;
  _nextConcentration = _concentration * decayFactor (minutesPerTimeStep);
  //TODO: if concentration is negligible, return true to clean up this object
  return false;
}
