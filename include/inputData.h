/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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


//parse the xml scenario file
//

#ifndef INPUTDATABOINC_H
#define INPUTDATABOINC_H

#include "Global.h"
#include "util/BoincWrapper.h"
#include <scenario.hxx>
#include <string>
#include <bitset>

namespace OM {

namespace Params {
  enum Params {
    /// @b Used in NoVectorControl
    //@{
    NEG_LOG_ONE_MINUS_SINF = 1,
    E_STAR = 2,
    SIMM = 3,
    X_STAR_P = 4,
    GAMMA_P = 5,
    //@}
    SIGMA_I_SQ = 6,			///< Used in WithinHostModel
    /// @b Used in Infection
    //@{
    CUMULATIVE_Y_STAR = 7,
    CUMULATIVE_H_STAR = 8,
    NEG_LOG_ONE_MINUS_ALPHA_M = 9,
    DECAY_M = 10,
    SIGMA0_SQ = 11,
    X_NU_STAR = 12,
    //@}
    /// @b Used in PathogenesisModel
    //@{
    Y_STAR_SQ = 13,
    ALPHA = 14,
    //@}
    DENSITY_BIAS_NON_GARKI = 15,	///< Used in WithinHostModel
    BASELINE_AVAILABILITY_SHAPE = 16,	///< Used in InfectionIncidenceModel
    LOG_ODDS_RATIO_CF_COMMUNITY = 17,	///< Used in CaseManagementModel
    INDIRECT_RISK_COFACTOR = 18,	///< Used in PathogenesisModel
    NON_MALARIA_INFANT_MORTALITY = 19,	///< Used in Summary
    DENSITY_BIAS_GARKI = 20,		///< Used in WithinHostModel
    SEVERE_MALARIA_THRESHHOLD = 21,	///< Used in PathogenesisModel
    IMMUNITY_PENALTY = 22,		///< Used in WithinHostModel
    IMMUNE_EFFECTOR_DECAY = 23,		///< Used in WithinHostModel
    /// @b Used in PathogenesisModel
    //@{
    COMORBIDITY_INTERCEPT = 24,
    Y_STAR_HALF_LIFE = 25,
    Y_STAR_1 = 26,
    //@}
    ASEXUAL_IMMUNITY_DECAY = 27,	///< Used in WithinHostModel
    /// @b Used in PathogenesisModel
    //@{
    Y_STAR_0 = 28,
    
    CRITICAL_AGE_FOR_COMORBIDITY = 30,
    MUELLER_RATE_MULTIPLIER = 31,
    MUELLER_DENSITY_EXPONENT = 32,
    //@}
    MAX
  };
}

/** Used to describe which interventions are in use. */
namespace Interventions {
    enum Flags {
	CHANGE_HS,
	CHANGE_EIR,
	VACCINE,	// any vaccine
	MDA,
	IPTI,
	ITN,
	IRS,
	VEC_AVAIL,
	LARVICIDING,
	SIZE
    };
}
    
    
    class InputDataType {
    public:
	InputDataType () : scenario(NULL) {}
	
	/** @brief Reads the document in the xmlFile
	* 
	* Throws on failure. */
	util::Checksum createDocument(std::string);

	/**
	* Some elements in memory have been created. This function deletes the object in memory
	*/
	void cleanDocument();
	
	/// Get the base scenario element
	const scnXml::Scenario& getScenario();
	
	/// Get the Monitoring xml object
	const scnXml::Monitoring& getMonitoring();
	/// Get the Interventions xml object
	const scnXml::Interventions& getInterventions();
	/// Get the EntoData xml object
	const scnXml::EntoData& getEntoData();
	/// Get the Demography xml object
	const scnXml::Demography& getDemography();
	/// Get the EventScheduler xml object
	const scnXml::EventScheduler& getEventScheduler();
	/// Get the HealthSystem xml object
	const scnXml::HealthSystem& getHealthSystem();

	/** Get a mutable version of scenario element.
	*
	* This is the only entry point for changing the scenario document.
	* 
	* You should set "documentChanged = true;" if you want your changes saved. */
	scnXml::Scenario& getMutableScenario();

	/** Change the element returned by getHealthSystem() to newHS.
	*
	* This doesn't change the scenario document, but just a local pointer, so it's
	* safe to use when editing and saving the document. */
	void changeHealthSystem (const scnXml::HealthSystem* newHS);

	/// Get the intervention from interventions->timed with time time.
	/// @returns NULL if not available
	const scnXml::Intervention* getInterventionByTime(int time);

	/// Returns and enum representing which interventions are active.
	const bitset<Interventions::SIZE> getActiveInterventions ();


	/// Get a parameter from the parameter list. i should be less than Params::MAX.
	double getParameter (size_t i);

	/// Set true if the xml document has been changed and should be saved.
	bool documentChanged;

	// -----  Other parameter-getters (old functions)  -----

	// For WHM:
	double get_detectionlimit(); 
	int get_analysis_no(); 

	// For global / simulation:
	double get_maximum_ageyrs();
	int get_latentp(); 
	int get_interval(); 

	// For population:
	int get_populationsize(); 

	// For transmission:
	int get_mode(); 
	double get_growthrate(); 

	// For GSL:
	int getISeed(); 
	
    private:
	void initParameterValues ();
	void initTimedInterventions ();
	
	/// Sometimes used to save changes to the xml.
	std::string xmlFileName;
	
	/** @brief The xml data structure. */
	scnXml::Scenario* scenario;
	const scnXml::Monitoring * monitoring;
	const scnXml::Interventions * interventions;
	const scnXml::EntoData * entoData; // May be replaced by a changeEIR intervention
	const scnXml::Demography * demography;
	const scnXml::HealthSystem * healthSystem; // May be replaced by a changeHS intervention or not present
	const scnXml::EventScheduler * eventScheduler; // Optional (may be NULL)
	const scnXml::Parameters * parameters;
	
	// Initialized (derived) values:
	double parameterValues[Params::MAX];
	std::map<int, const scnXml::Intervention*> timedInterventions;
	bitset<Interventions::SIZE> activeInterventions;
    };
    /// InputData entry point. Most set up done before sim start, but changeHealthSystem() and documentChanged members allow changes which aren't checkpointed. TODO
    extern InputDataType InputData;
}
#endif
