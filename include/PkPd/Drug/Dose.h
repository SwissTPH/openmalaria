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

/** A simple class to hold dose info. */
class Dose {
public:
  /** Create a new dose. */
  Dose () : time(0.0), mg(0.0) {}
  /// ditto
  Dose (double time, double mg) {
    this->time = time;
    this->mg = mg;
  }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
	time & stream;
	mg & stream;
    }
  
  /// time (days)
  double time;
  /// dose in mg
  double mg;
};
} }
#endif
