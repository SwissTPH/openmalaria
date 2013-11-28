/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2013 Swiss Tropical and Public Health Institute 
 * Copyright (C) 2005-2013 Liverpool School Of Tropical Medicine
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

#include "WithinHost/ProphylacticActionWithinHost.h"

#include "util/random.h"
#include <boost/concept_check.hpp>

namespace OM { namespace WithinHost {
    using namespace util;

void ProphylacticActionWithinHost::addProphylacticEffects(
    const vector<double>& pClearanceByTime)
{
    Pending_t::iterator pend_it = pendingClearanceProbabilities.begin();
    bool using_it = true;
    for( vector<double>::const_iterator it = pClearanceByTime.begin();
        it != pClearanceByTime.end(); ++it )
    {
        if( using_it && pend_it != pendingClearanceProbabilities.end() ){
            double p = pend_it->first;
            double n = pend_it->second;
            p = (n * p + 1 * *it) / (n + 1);
            pend_it->first = p;
            pend_it->second = n + 1;
            ++pend_it;
        }else{
            using_it = false;   // stop using the iterator (can we do so after appending elements to the list?)
            pendingClearanceProbabilities.push_back( pair<double,uint32_t>( *it, 1 ) );
        }
    }
}

void ProphylacticActionWithinHost::drugAction(){
    if( pendingClearanceProbabilities.empty() ) return;
    double pClearance = pendingClearanceProbabilities.front().first;
    if( util::random::bernoulli( pClearance ) )
        clearAllInfections();
    pendingClearanceProbabilities.pop_front();
}

// -----  Data checkpointing  -----

void ProphylacticActionWithinHost::checkpoint (istream& stream) {
    DescriptiveWithinHostModel::checkpoint (stream);
    size_t len;
    len & stream;
    validateListSize(len);
    for( size_t i =0; i < len; ++i ){
        double p;
        p & stream;
        uint32_t n;
        n & stream;
        pendingClearanceProbabilities.push_back( make_pair( p, n ) );
    }
}
void ProphylacticActionWithinHost::checkpoint (ostream& stream) {
    DescriptiveWithinHostModel::checkpoint (stream);
    pendingClearanceProbabilities.size() & stream;
    for( Pending_t::const_iterator it = pendingClearanceProbabilities.begin();
        it != pendingClearanceProbabilities.end(); ++it )
    {
        it->first & stream;
        it->second & stream;
    }
}

} }