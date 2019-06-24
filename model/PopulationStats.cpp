/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

#include "PopulationStats.h"

namespace OM {
    boost::int64_t PopulationStats::totalInfections =0;
    boost::int64_t PopulationStats::allowedInfections =0;
    boost::int64_t PopulationStats::humanUpdateCalls =0;
    boost::int64_t PopulationStats::humanUpdates =0;
    
    void PopulationStats::print() {
	long double x = 100.0 * (totalInfections - allowedInfections) / totalInfections;
	cerr
	    << "Total/allowed infections: "
	    <<totalInfections<<"/"<<allowedInfections
	    <<"\t("<<x<<"% prevented)"
	    <<endl
	;
	
	x = 100.0 * (humanUpdateCalls - humanUpdates) / humanUpdateCalls;
	cerr
	    << "Human updates/total calls: "
	    <<humanUpdates<<"/"<<humanUpdateCalls
	    <<"\t("<<x<<"% skipped)"
	    <<endl
	;
    }
    
    void PopulationStats::staticCheckpoint (istream& stream){
	totalInfections & stream;
	allowedInfections & stream;
	humanUpdateCalls & stream;
	humanUpdates & stream;
    }
    void PopulationStats::staticCheckpoint (ostream& stream){
	totalInfections & stream;
	allowedInfections & stream;
	humanUpdateCalls & stream;
	humanUpdates & stream;
    }
    
}
