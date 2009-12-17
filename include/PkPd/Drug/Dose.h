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

#ifndef Hmod_Dose
#define Hmod_Dose

#include "Global.h"
#include "PkPd/Proteome.h"

using namespace std;

namespace OM { namespace PkPd {

/** A simple class to hold dose info.
 *
 * TODO: This was created as a place-holder. Replace x,y with whatever data is
 * needed, and remove what isn't. */
class Dose {
public:
  /** Create a new dose. */
  Dose (double x, double y) {
    this->x = x;
    this->y = y;
  }
  /** Load from a checkpoint. */
  Dose (istream& in) {
    in >> x;
    in >> y;
  }
  
  /** Write a checkpoint. */
  void write (ostream& out) const {
    out << x << endl;
    out << y << endl;
  }
  
  /// Some type of data is wanted... (concentration at start of next timestep, and integral of concentration for this timestep?)
  double x,y;
};
} }
#endif
