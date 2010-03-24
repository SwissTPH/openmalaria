/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#ifndef Hmod_PopulationAgeStructure
#define Hmod_PopulationAgeStructure

#include <vector>

/*  TODO:
    Provide saving and loading functionality for cumAgeProp, to allow building a
    version of the program without GSL (useful if random-number generators are
    also re-written).

1.  Needs to be stored somewhere; probably <demography> in the XML document
    is a good place (InputData already can save changes to the scenario document).
2.  If compiled with GSL, calculate cumAgeProp (as currently done) and save/
    validate. If not, use cumAgeProp from the XML (or complain if not present).
3.  OM_HAVE_GSL can become a macro added by CMake.
******************************************************************************/
#define OM_HAVE_GSL

namespace OM
{
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
#ifdef OM_HAVE_GSL
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
#endif	// defined OM_HAVE_GSL
	
	//BEGIN static parameters set by init() or calcCumAgeProp()
        //! max lifespan in intervals
	static int maxTimestepsPerLife;
	
	/** Target cumulative percentage of population by age, from oldest age to youngest.
	*
	* cumAgeProp[_maxTimestepsPerLife+1-i] gives the proportion of people aged i timesteps or older.
	*/
	static std::vector<double> cumAgeProp;
	//END
    };
}

#endif
