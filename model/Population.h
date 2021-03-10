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
    
    /// Type of population list. Store pointers to humans only to avoid copy operations.
    typedef vector<Host::Human> HumanPop;
    /// Iterator type of population
    typedef HumanPop::iterator Iter;
    /// Const iterator type of population
    typedef HumanPop::const_iterator ConstIter;
    /// Const reverse iterator type of population
    typedef HumanPop::const_reverse_iterator ConstReverseIter;
    
    /** @brief Access the population list, as a whole or with iterators. */
    //@{
    // non-const versions are needed to do things like add infections
    inline Iter begin() { return population.begin(); }
    inline Iter end() { return population.end(); }
    inline ConstIter cbegin() const{ return population.cbegin(); }
    inline ConstIter cend() const{ return population.cend(); }
    inline ConstReverseIter crbegin() const{ return population.crbegin(); }
    inline ConstReverseIter crend() const{ return population.crend(); }
    // Pair of iterators (begin, end)
    inline std::pair<Iter, Iter> range() {
        return std::make_pair(population.begin(), population.end());
    }
    // Pair of const iterators (cbegin, cend)
    inline std::pair<ConstIter, ConstIter> crange() const {
        return std::make_pair(population.cbegin(), population.cend());
    }
    /** Return the number of humans. */
    inline size_t size() const {
        return populationSize;
    }
    inline HumanPop &getHumans() {
        return population;
    }
    inline const HumanPop &getHumans() const {
        return population;
    }
    //@}

private:
    /// Delegate to print the number of hosts
    void ctsHosts (ostream& stream);
    /// Delegate to print cumulative numbers of hosts under various age limits
    void ctsHostDemography (ostream& stream);
    /// Delegate to print the number of births since last count
    void ctsRecentBirths (ostream& stream);
    /// Delegate to print the number of patent hosts
    void ctsPatentHosts (ostream& stream);
    /// Delegate to print immunity's cumulativeh parameter
    void ctsImmunityh (ostream& stream);
    /// Delegate to print immunity's cumulativeY parameter (mean across population)
    void ctsImmunityY (ostream& stream);
    /// Delegate to print immunity's cumulativeY parameter (median across population)
    void ctsMedianImmunityY (ostream& stream);
    /// Delegate to print the mean age-based availability reduction of each human relative to an adult
    void ctsMeanAgeAvailEffect (ostream& stream);
    void ctsITNCoverage (ostream& stream);
    void ctsIRSCoverage (ostream& stream);
    void ctsGVICoverage (ostream& stream);
    /// Delegate to print the mean hole index of all bed nets
//     void ctsNetHoleIndex (ostream& stream);
    

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
    HumanPop population;
    
    friend class AnophelesModelSuite;
};

}
#endif
