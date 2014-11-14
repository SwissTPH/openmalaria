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

#ifndef H_OM_mon_cpp
// This is the only way we need to use this header:
#error "Only include from mon.cpp"
#endif

#include "mon/reporting.h"      // for Measure enum
#include <map>
#include <string>

namespace OM {
namespace mon {

// Describes each "measure" to be output
struct OutMeasure{
    int outId;  // number used in output to identify this measure/aggregation
    Measure m;  // internal measure (e.g. MHR_HOSTS) this comes from
    bool isDouble;  // false: type is int; true: type is double
    bool byAge; // segregate by age
    bool byCohort;      // segregate by cohort
    
    // Convenience constructors:
    OutMeasure() : outId(-1), m(M_NUM), isDouble(false), byAge(false), byCohort(false) {}
    OutMeasure( int outId, Measure m, bool isDouble, bool byAge, bool byCohort ) :
        outId(outId), m(m), isDouble(isDouble), byAge(byAge), byCohort(byCohort)
    {}
    // Something with reports segregated by human age and cohort membership
    static OutMeasure humanAC( int outId, Measure m, bool isDouble ){
        return OutMeasure( outId, m, isDouble, true, true );
    }
};

// These are all output measures set by a name in the XML
// Example: nHosts
typedef std::map<std::string,OutMeasure> NamedMeasureMapT;
NamedMeasureMapT namedOutMeasures;

// This method defines output measures accepted by name in the XML (e.g.
// "nHost") and their numeric output identifier (i.e. measure column of
// outputs), type of output (integer or floating point), aggregation, and the
// corresponding internal measure code.
void defineOutMeasures(){
    /// Total number of humans
    namedOutMeasures["nHost"] = OutMeasure::humanAC( 0, MHR_HOSTS, false );
    /** The number of human hosts with an infection (patent or not) at the time
     * the survey is taken. */
    namedOutMeasures["nInfect"] = OutMeasure::humanAC( 1, MHR_INFECTED_HOSTS, false );
    /** Expected number of infected hosts
     * 
     * This is the sum of the probabilities, across all time steps since the
     * last survey, of each host becoming infected on that time step. */
    namedOutMeasures["nExpectd"] = OutMeasure::humanAC( 2, MHD_EXPECTED_INFECTED, true );
    /** The number of human hosts whose total (blood-stage) parasite density is
     * above the detection threshold */
    namedOutMeasures["nPatent"] = OutMeasure::humanAC( 3, MHR_PATENT_HOSTS, false );
    /** The total number of infections in the population: includes both blood
     * and liver stages. Vivax: this is the number of broods. */
    namedOutMeasures["totalInfs"] = OutMeasure::humanAC( 6, MHR_INFECTIONS, false );
    /** The sum of all detectable infections (where blood stage parasite
     * density is above the detection limit) across all human hosts.
     * Vivax: the number of broods with an active blood stage. */
    namedOutMeasures["totalPatentInf"] = OutMeasure::humanAC( 8, MHR_PATENT_INFECTIONS, false );
    
    /** Sum (across hosts) of the natural logarithm of the parasite density of
     * hosts with detectable parasite density (patent according to the
     * monitoring diagnostic). */
    namedOutMeasures["sumlogDens"] = OutMeasure::humanAC( 5, MHD_LOG_DENSITY, true );
    /// Sum of log(1 + p) where p is the pyrogenic threshold
    namedOutMeasures["sumLogPyrogenThres"] = OutMeasure::humanAC( 4, MHD_LOG_PYROGENIC_THRESHOLD, true );
    /// Sum of the pyrogenic threshold
    namedOutMeasures["sumPyrogenThresh"] = OutMeasure::humanAC( 10, MHD_PYROGENIC_THRESHOLD, true );
}

}
}
