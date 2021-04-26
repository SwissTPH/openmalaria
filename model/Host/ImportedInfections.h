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

#ifndef OM_IMPORTED_INFECTIONS
#define OM_IMPORTED_INFECTIONS

#include "Global.h"
#include "schema/interventions.h"
#include "util/errors.h"

namespace OM {
    class Population;

namespace Host {

    class ImportedInfections {
    public:
        ImportedInfections() : period(sim::zero()), lastIndex(0) {}
        
        /** Initialise, passing intervention description
         * 
         * @param iiElt The scenario element with description of rates
         * @returns True if any infections are imported
         */
        void init( const scnXml::ImportedInfections& iiElt );
        
        /** Import this time-step's imported infections, according to initialised rates
         * 
         *  The probability of an host to import an infection is calculated
         *  from the importedInfectionsPerThousandHosts. The bernoulli distribution
         *  is then used to predict if an human has imported the infection in the
         *  population or not. A maximum of one infection can be imported per
         *  person.
         * 
         * @param pop The Population class encapsulating all humans */
        void import( Population& pop );
        
        /// Checkpointing
        template<class S>
        void operator& (S& stream) {
            using namespace OM::util::checkpoint;
            // period and rate are set from XML and not changed
            lastIndex & stream;
        }
        
    private:
        SimTime period;
        uint32_t lastIndex;
        struct Rate {
            Rate( SimTime t, double v ): time(t), value(v) {}
            SimTime time;
            double value;
            inline bool operator< (const Rate& that) const{
                return time < that.time;
            }
        };
        vector<Rate> rate;
    };
} }

#endif
