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

#ifndef Hmod_dose
#define Hmod_dose

//#include <string>
//#include <deque>
//#include <map>
//#include <vector>
#include "Global.h"
#include "proteome.h"

using namespace std;

/** A simple class to hold dose info. */
class Dose {
public:
  /** Create a new dose. */
  Dose (double x, double y) {
    this->x = x;
    this->y = y;
  }
  /** Load from a checkpoint. */
  Dose (istream& in);
  /** Write a checkpoint. */
  void write (ostream& out) const;
  
  /// Some type of data is wanted... (concentration at start of next timestep, and integral of concentration for this timestep?)
  double x,y;
};

#endif
