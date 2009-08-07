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


//parse the xml scenario file
//

#ifndef INPUTDATABOINC_H
#define INPUTDATABOINC_H

#include "Constant.h"
#include <scenario.hxx>
#include <string>
class Anopheles;


/** @brief Reads the document in the xmlFile
 * 
 * Throws on failure. */
void createDocument(std::string);

/**
* Some elements in memory have been created. This function deletes the object in memory
*/
void cleanDocument();

/// Get the Monitoring xml object
const scnXml::Monitoring& getMonitoring();
/// Get the Interventions xml object
const scnXml::Interventions& getInterventions();
/// Get the EntoData xml object
const scnXml::EntoData& getEntoData();
/// Get the Demography xml object
const scnXml::Demography& getDemography();
/// Get the CaseManagements xml object
const scnXml::CaseManagements* getCaseManagements();
/// Get the HealthSystem xml object
const scnXml::HealthSystem& getHealthSystem();

// Change xml data for certain interventions
void changeHealthSystem (const scnXml::HealthSystem*);

/// Get the intervention from interventions->timed with time time.
/// @returns NULL if not available
const scnXml::Intervention* getInterventionByTime(int time);

/// Get a parameter from the parameter list. i should be less than Params::MAX.
double getParameter (size_t i);

/// Set true if the xml document has been changed and should be saved.
extern bool documentChanged;

// -----  Other parameter-getters (old functions)  -----

// For WHM:
double get_detectionlimit(); 
int get_analysis_no(); 

// For summary:
int get_number_of_surveys(); 
int get_time_of_survey(int ); 
int get_summary_option();
int get_assim_mode(); 
double get_lowerbound(); 

// For global / simulation:
int get_model_version();
double get_maximum_ageyrs();
int get_simulation_duration(); 
int get_latentp(); 
int get_interval(); 

// For population:
int get_wu_id();
int get_populationsize(); 

// For transmission:
int get_mode(); 
double get_growthrate(); 

// For GSL:
int getISeed(); 

#endif
