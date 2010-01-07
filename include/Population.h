/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include "Host/Human.h"
#include "Transmission/TransmissionModel.h"
#include "inputData.h"

#include <list>
#include <fstream>
#include <vector>

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

    /** Checks for time-based interventions and implements them
     *
     * \param time Current time (in tsteps) */
    void implementIntervention (int time);

private:
    //! Creates initializes and add to the population list a new uninfected human
    /*!
       \param dob date of birth (usually current time)
    */
    void newHuman (int dob);
    
    /** Generic function to activate some intervention on all humans within the
     * age range and passing the compliance test given by mass.
     *
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention A member-function pointer to a "void func ()" function
     * within human which activates the intervention. */
    void massIntervention (const scnXml::Mass& mass, void (Host::Human::*intervention) ());

    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);


    //! Size of the human population
    int populationSize;

public:
    //! TransmissionModel model
    Transmission::TransmissionModel* _transmissionModel;

private:
    /** The simulated human population
     *
     * The list of all humans, ordered from oldest to youngest. */
    std::list<Host::Human> population;

    /// Iterator type of population
    typedef std::list<Host::Human>::iterator HumanIter;

    friend class VectorAnophelesSuite;


    /** Encapsulates code just setting up the age structure (i.e.  cumAgeProp). */
    class AgeStructure
    {
    public:
	/** Set up cumAgeProp from XML data. */
	static void init ();
	
	/** Return maximum individual lifetime in intervals that AgeStructure can handle. */
	static inline int getMaxTimestepsPerLife () {
	    return maxTimestepsPerLife;
	}
	
	/** Return the expected population size of individuals aged ageTSteps or
	* older, based on a total population size of targetPop. */
	static int targetCumPop (int ageTSteps, int targetPop);
	
    private:
        /*! Estimates demography parameters to define a smooth curve for the target
        population age-distribution (age in years) */
        static void estimateRemovalRates();
	
	static double minimizeCalc_rss(double* par1, double* par2);
	
        /** For input values for alpha1 and mu1, the fit to field data (residualSS)
        * is calculated and returned function called iteratively by
        * estimateRemovalRates. */
        static double setDemoParameters (double param1, double param2);

        /** Takes the best-fitting demography parameters estimated by
        * estimateRemovalRates and calculates the age structure (cumAgeProp). */
	static void calcCumAgeProp ();

        /** This is the maximum age of an individual that the simulation program can
        * handle. Max age for a scenario is given in the  the xml file. */
        static const int maxLifetimeDays = 32855;
        static const int ngroups = 20;
	
	//BEGIN static parameters set by init() or calcCumAgeProp()
        //! max lifespan in intervals
	static int maxTimestepsPerLife;
	
	/** Target cumulative percentage of population by age, from oldest age to youngest.
	*
	* cumAgeProp[_maxTimestepsPerLife+1-i] gives the proportion of people aged i timesteps or older.
	*/
	static vector<double> cumAgeProp;
	//END
	
	//BEGIN static parameters only used by estimateRemovalRates(), setDemoParameters() and calcCumAgeProp()
	/** The bounds for each age group and percentage of population in this age
        * group for the field data demography age groups.
        *
        * ageGroupBounds[i] is the lower-bound for group i, ageGroupBounds[i+1] is
        * the group's upper bound. ageGroupPercent[i] is the percentage of the
        * population in age group i.
        *
        * Set by estimateRemovalRates() and used internally (by
        * setDemoParameters()). */
        //@{
        static double ageGroupBounds[ngroups+1];
	static double ageGroupPercent[ngroups];
        //@}
        /** Demography variables used in estimating the smooth curve.
        *
        * Only used in setDemoParameters() calculations. */
        //@{
        static double M1[ngroups];
	static double M2[ngroups];
	static double M[ngroups];
	static double pred[ngroups];
        //@}
        /** Parameters defining smooth curve of target age-distribution.
        *
        * Set by estimateRemovalRates() (via setDemoParameters()) and used by
        * setupPyramid(). */
        //@{
	static double mu0;
	static double mu1;
	static double alpha0;
	static double alpha1;
	static double rho;
        //@}
	//END
    };
};

}
#endif
