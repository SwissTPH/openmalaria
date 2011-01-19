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
#include "util/AgeGroupInterpolation.h"

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
	/** Initialisation - reads fatality rates, etc.
	 * Both derived case management systems have their own init function;
	 * this is named "initCommon" to avoid confusion over which is called. */
	static void initCommon ();
        /// Free memory
        static void cleanupCommon ();
	
	/** Checkpointing: load state.
	 *
	 * If healthSystemSource != TimeStep::never, this calls changeHealthSystem to re-load
	 * a parameters from an intervention. */
	static void staticCheckpoint (istream& stream);
	/// Checkpointing: save state
	static void staticCheckpoint (ostream& stream);
	
	/** Serves both to set initial health-system data and to change
	 * following an intervention.
	 * 
	 * Gets the primary description of the health system when
	 * source == TimeStep::never, or the replacement description given by
	 * timed intervention at timestep source, then calls
	 * setHealthSystem from a derived class.
	 * 
	 * Also calls readCaseFatalityRatio with the new data. */
	static void changeHealthSystem (TimeStep source);
	
	/// Return the case-fatality-rate map (needed by EventScheduler)
// 	static inline const map<double,double>& getCaseFatalityRates (){
// 	    return caseFatalityRates;
// 	}
	
	/** Calculate the case fatality rate in the community as a function of
	 * the hospital case fatality rate. */
	static double getCommunityCaseFatalityRate(double caseFatalityRatio);
	
	/** Stepwise linear interpolation to get age-specific hospital case
	 * fatality rate from input data.
	 * 
	 * @param ageYears Age of person in years */
	static inline double caseFatality(double ageYears) {
            return caseFatalityRate->eval( ageYears );
        }
        
        /** Scale case fatality rates by a factor. */
        static inline void scaleCaseFatalityRate( double factor ){
            caseFatalityRate->scale( factor );
        }
	
	/** Get the probability of in-hospital sequale for a severe bout
	 * (step-wise constant; I don't think anything else makes sense since
	 * data is currently from two age-groups).
	 * 
	 * Currently we use the same values for outpatients. */
	static inline double pSequelaeInpatient(double ageYears) {
            return pSeqInpatient->eval(ageYears);
        }
	
    private:
	/** Gets the primary description of the health system when
	 * source == TimeStep::never, or the replacement description given by
	 * timed intervention at timestep source. */
	static const scnXml::HealthSystem& getHealthSystem ();
	
	/// Reads the CFR and sequelae data from healthSystem.
	static void readCommon(const scnXml::HealthSystem& healthSystem);
	
	/** Describes which health-system descriptor should be used, in order
	 * to load the correct one from a checkpoint (see getHealthSystem). */
	static TimeStep healthSystemSource;
	
	// An interpolation function from age-groups to CFR values.
        static util::AgeGroupInterpolation* caseFatalityRate;
	
	//log odds ratio of case-fatality in community compared to hospital
	static double _oddsRatioThreshold;
	
	static util::AgeGroupInterpolation* pSeqInpatient;
      };
} }
#endif
