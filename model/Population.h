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
#include <fstream>
#include <utility>  // pair

namespace scnXml{
    class Scenario;
}
namespace OM {
    class Parameters;

//! The simulated human population
class Population
{
public:
    /// Call static inits of sub-models
    static void init( const OM::Parameters& parameters, const scnXml::Scenario& scenario );

    /// Checkpointing for static data members
    static void staticCheckpoint (istream& stream);
    static void staticCheckpoint (ostream& stream); ///< ditto


    Population( size_t populationSize );
    
    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);
    
    /** Creates the initial population of Humans according to cumAgeProp. */
    void createInitialHumans();
    
    /** Initialisation run between initial one-lifespan run of simulation and
     * actual simulation. */
    void preMainSimInit ();

    //! Updates all individuals in the list for one time-step
    /*!  Also updates the population-level measures such as infectiousness, and
         the age-distribution by c outmigrating or creating new births if
         necessary */
    void update( Transmission::TransmissionModel& transmission, SimTime firstVecInitTS );

    //! Makes a survey
    void newSurvey();
    
    /// Flush anything pending report. Should only be called just before destruction.
    void flushReports();

    //! Size of the human population
    size_t populationSize;
    
    ///@brief Variables for continuous reporting
    //@{
    vector<double> ctsDemogAgeGroups;
    
    /// Births since last continuous output
    int recentBirths;
    //@}
    
    /** The simulated human population
     *
     * The list of all humans, ordered from oldest to youngest. */
    vector<Host::Human> humans;
};

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
/// Delegate to print the mean hole index of all bed nets
//     void ctsNetHoleIndex (ostream& stream);

}
#endif
