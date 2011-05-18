/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include "Transmission/TransmissionModel.h"

#include <list>
#include <fstream>

namespace OM
{

//! The simulated human population
class Population
{
public:
    /// Call static inits of sub-models
    static void init();

    /// Calls static clear on sub-models to free memory
    static void clear();

    /// Checkpointing for static data members
    static void staticCheckpoint (istream& stream);
    static void staticCheckpoint (ostream& stream); ///< ditto


    Population();
    //! Clears human collection.
    ~Population();

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        populationSize & stream;
	recentBirths & stream;
        (*_transmissionModel) & stream;
	
        checkpoint (stream);
    }
    
    
    /** Creates the initial population of Humans according to cumAgeProp.
    *
    * Also runs some transmission model initialisation (which needs to happen
    * after a population has been created). */
    void createInitialHumans ();
    
    /** Initialisation run between initial one-lifespan run of simulation and
     * actual simulation. */
    void preMainSimInit ();

    //! Updates all individuals in the list for one time-step
    /*!  Also updates the population-level measures such as infectiousness, and
         the age-distribution by c outmigrating or creating new births if
         necessary */
    void update1();

    //! Makes a survey
    void newSurvey();
    
    /// Flush anything pending report. Should only be called just before destruction.
    void flushReports();
    
    /** Checks for time-based interventions and implements them
     *
     * \param time Current time (in tsteps) */
    void implementIntervention (TimeStep time);
    
    /// Type of population list
    typedef list<Host::Human> HumanPop;
    /// Iterator type of population
    typedef HumanPop::iterator HumanIter;
    
    /** Return the list of humans. */
    inline HumanPop& getList() {
        return population;
    }

private:
    //! Creates initializes and add to the population list a new uninfected human
    /*!
       \param dob date of birth (usually current time)
    */
    void newHuman (TimeStep dob);
    
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
    
    /** This function sets the imported infections in a population.
     *  The probability of an host to import an infection is calculated
     *  from the importedInfectionsPerThousandHosts. The bernoulli distribution
     *  is then used to predict if an human has imported the infection in the
     *  population or not. A maximum of one infection can be imported per
     * intervention.
     * 
     * @param rate Number of infections imported per 1000 humans at this time
     * point.  */
    void importedInfections(double rate);

    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);


    //! Size of the human population
    int populationSize;
    
    ///@brief Variables for continuous reporting
    //@{
    vector<double> ctsDemogAgeGroups;
    
    /// Births since last continuous output
    int recentBirths;
    //@}
public:
    //! TransmissionModel model
    Transmission::TransmissionModel* _transmissionModel;
    
private:
    /** The simulated human population
     *
     * The list of all humans, ordered from oldest to youngest. */
    HumanPop population;

    friend class VectorAnophelesSuite;
};

}
#endif
