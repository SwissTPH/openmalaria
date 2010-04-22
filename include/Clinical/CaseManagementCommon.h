/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef Hmod_CaseManagementCommon
#define Hmod_CaseManagementCommon

#include "Global.h"
#include <map>

namespace scnXml {
    class HealthSystem;
}

namespace OM { namespace Clinical {
    /** Code shared between case management systems, to:
     *
     * 1. deal with change-of-health-system interventions
     * 2. handle case-fatality-rate data
     * 
     * For convenience with handling change-of-health-system interventions,
     * the case management model holds CFR data, not the clinical model.
     *************************************************************************/
    class CaseManagementCommon {
    public:
	/** Checkpointing: load state.
	 *
	 * If healthSystemSource != -1, this calls changeHealthSystem to re-load
	 * a parameters from an intervention. */
	static void staticCheckpoint (istream& stream);
	/// Checkpointing: save state
	static void staticCheckpoint (ostream& stream);
	
	/** Serves both to set initial health-system data and to change
	 * following an intervention.
	 * 
	 * Gets the primary description of the health system when
	 * source == -1, or the replacement description given by
	 * timed intervention at timestep source, then calls
	 * setHealthSystem from a derived class.
	 * 
	 * Also calls readCaseFatalityRatio with the new data. */
	static void changeHealthSystem (int source);
	
	/// Return the case-fatality-rate map (needed by EventScheduler)
	static inline const map<double,double>& getCaseFatalityRates (){
	    return caseFatalityRates;
	}
	
    protected:
	/** Stepwise linear interpolation to get age-specific hospital case
	 * fatality rate from input data.
	 * 
	 * @param ageYears Age of person in years */
	static double caseFatality(double ageYears);
	
    private:
	/** Gets the primary description of the health system when
	 * source == -1, or the replacement description given by
	 * timed intervention at timestep source. */
	static const scnXml::HealthSystem& getHealthSystem ();
	
	/// Reads in the Case Fatality percentages from the XML.
	static void readCaseFatalityRatio(const scnXml::HealthSystem& healthSystem);
	
	/** Describes which health-system descriptor should be used, in order
	 * to load the correct one from a checkpoint (see getHealthSystem). */
	static int healthSystemSource;
	
	// A map from age-group upper-bounds to CFR values. The first lower-
	// bound we assume to be less than or equal to any input value, and we
	// assume no input can be greater than the last upper-bound (which
	// should be inf).
	static map<double,double> caseFatalityRates;
    };
} }
#endif
