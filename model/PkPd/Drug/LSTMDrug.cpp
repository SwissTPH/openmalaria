/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#include "PkPd/Drug/LSTMDrug.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <algorithm>

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrug::LSTMDrug()
{}

// Create our list of doses. Optimise for the case where we only have 1 or less per day, but be able to handle more.
// If two oral doses coincide, they can be combined but doing so is not essential.
// End of an IV dose can be represented by an oral dose of 0 qty potentially.
// IV doses must be split over the end of the day and if an oral dose occurs in the middle.
// Overlapping IV doses are not supported.

void LSTMDrug::_medicate (double time, double qty, double volDist) {
    double conc = qty / volDist;        // mg / l
    // multimap insertion: is ordered
    DoseMap::iterator lastInserted =
    doses.insert (doses.end(), make_pair (time, DoseParams( conc, 0 )));
    check_split_IV( lastInserted );
}

void LSTMDrug::medicateIV (double time, double duration, double qty ) {
    assert( duration > 0.0 );
    
    double infusRate = qty / duration;	// mg/day
    DoseMap::iterator lastInserted =
    doses.insert(doses.end(), make_pair( time, DoseParams( infusRate, duration ) ) );    
    check_split_IV( lastInserted );
    
    doses.insert( doses.end(), make_pair( time+duration, DoseParams() ) );
}

void LSTMDrug::check_split_IV( DoseMap::iterator lastInserted ){
    for( DoseMap::iterator it = doses.begin(); it != doses.end(); ++it ){
        if( it == lastInserted )
            continue;
        
        if( it->first == lastInserted->first ){
            if( it->second.duration == 0.0 && lastInserted->second.duration == 0.0 ){
                it->second.qty += lastInserted->second.qty;
                doses.erase( lastInserted );
                return;
            } else if( it->second.duration == 0.0 ){
                // oral followed by IV; no problem
            } else if( lastInserted->second.duration == 0.0 ){
                // IV followed by oral; needs to be analysed in other order
                swap( it->second, lastInserted->second );
            } else {
                throw util::xml_scenario_error( "IV,IV medications overlap — not supported!" );
            }
        } else if( it->first < lastInserted->first ){
            if( it->first + it->second.duration > lastInserted->first ){
                // TODO: if one is an oral dose, could split IV
                throw util::xml_scenario_error( "IV/oral medications overlap — not supported!" );
            }
        } else if( it->first > lastInserted->first ){
            if( lastInserted->first + lastInserted->second.duration > it->first ){
                // TODO: if one is an oral dose, could split IV
                throw util::xml_scenario_error( "IV/oral medications overlap — not supported!" );
            }
        }
    }
    
    // Make sure one dose is evaluated separately on each day
    double dayEnd = 1.0;
    while( lastInserted->first + lastInserted->second.duration > dayEnd ){
        double remaining_duration = lastInserted->first + lastInserted->second.duration - dayEnd;
        lastInserted->second.duration = dayEnd - lastInserted->first;
        lastInserted = doses.insert( doses.end(), make_pair( dayEnd, DoseParams( lastInserted->second.qty, remaining_duration ) ) );
        dayEnd += 1.0;
    }
}

}

namespace util { namespace checkpoint {
    
void operator& (multimap<double,PkPd::DoseParams> x, ostream& stream) {
    x.size() & stream;
    for (multimap<double,PkPd::DoseParams>::iterator pos = x.begin (); pos != x.end() ; ++pos) {
        pos->first & stream;
        pos->second & stream;
    }
}
void operator& (multimap<double,PkPd::DoseParams>& x, istream& stream) {
    size_t l;
    l & stream;
    validateListSize (l);
    x.clear ();
    multimap<double,PkPd::DoseParams>::iterator pos = x.begin ();
    for (size_t i = 0; i < l; ++i) {
        double s;
        PkPd::DoseParams t;
        s & stream;
        t & stream;
        pos = x.insert (pos, make_pair (s,t));
    }
    assert( x.size() == l );
}

} }

}