/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
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

#ifndef Hmod_Population
#define Hmod_Population

#include "Global.h"
#include "PopulationAgeStructure.h"
#include "Host/Human.h"

#include <vector>

namespace OM {

//! The simulated human population
class Population
{
public:
    Population(size_t size);

    /** Creates the initial population of Humans */
    void createInitialHumans();

    /** Remove dead humans, outmigrate others and introduce babies to
     * while keeping the population size and demography distirbution unchanged. */
    void regularize();

    /** Returns the size of the human population */
    size_t getSize() const;

    /** Returns the number of recent births */
    int getRecentBirths() const;

    /** Reset the number of recent births to 0 */
    void resetRecentBirths();

    /** Checkpointing methods */
    void checkpoint(istream& stream);
    void checkpoint(ostream& stream);

    /** The simulated human population */
    vector<Host::Human> humans;

private:
    /** Size of the human population */
    size_t size = 0;
    
    /** Births since last continuous output */
    int recentBirths = 0;
};

void registerContinousPopulationCallbacks();

/// Delegate to print the number of hosts
void ctsHosts (Population &population, ostream& stream);
/// Delegate to print cumulative numbers of hosts under various age limits
void ctsHostDemography (Population &population, ostream& stream);
/// Delegate to print the number of births since last count
void ctsRecentBirths (Population &population, ostream& stream);
/// Delegate to print the number of patent hosts
void ctsPatentHosts (Population &population, ostream& stream);
/// Delegate to print immunity's cumulativeh parameter
void ctsImmunityh (Population &population, ostream& stream);
/// Delegate to print immunity's cumulativeY parameter (mean across population)
void ctsImmunityY (Population &population, ostream& stream);
/// Delegate to print immunity's cumulativeY parameter (median across population)
void ctsMedianImmunityY (Population &population, ostream& stream);
/// Delegate to print the mean age-based availability reduction of each human relative to an adult
void ctsMeanAgeAvailEffect (Population &population, ostream& stream);
void ctsITNCoverage (Population &population, ostream& stream);
void ctsIRSCoverage (Population &population, ostream& stream);
void ctsGVICoverage (Population &population, ostream& stream);

}
#endif
