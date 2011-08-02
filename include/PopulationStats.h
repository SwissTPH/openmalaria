/*
  This file is part of OpenMalaria.
 
  Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

/* This module is concerned with collecting population-level statistics for
 command-line output (i.e. principally debugging stuff). */

#ifndef Hmod_PopStats
#define Hmod_PopStats

#include "Global.h"

namespace OM {
    class PopulationStats {
    public:
	// Print collected information at end of simulation
	static void print ();
	
	/// Checkpointing for static data members
	static void staticCheckpoint (istream& stream);
	static void staticCheckpoint (ostream& stream); ///< ditto
	
	static boost::int64_t totalInfections;
	static boost::int64_t allowedInfections;
	
	static boost::int64_t humanUpdateCalls;
	static boost::int64_t humanUpdates;
    };
}

#endif

