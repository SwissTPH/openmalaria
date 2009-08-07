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


#include "inputData.h"
#include "scenario.hxx"
#include "global.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
using namespace std;
// This namespace should make it obvious when types come from scenario.hxx and
// when they come from model code. It is not explicitly used here (prepending
// scnXml::), but I recommend code within the model explicitly do this.
using namespace scnXml;

/// Current schema version.
const int SCHEMA_VERSION = 7;
/** Oldest which current code is potentially compatible with
 * (provided the scenario.xml file references this version and doesn't use
 * members changed in newer versions). */
const int OLDEST_COMPATIBLE = 5;
/// Sometimes used to save changes to the xml.
std::string xmlFileName;
/// Set true if the xml document is changed and should be saved
bool documentChanged = false;

/** @brief The xml data structure. */
const Scenario* scenario = NULL;
const Monitoring * monitoring;
const Interventions * interventions;
const EntoData * entoData;	// May be replaced by a changeEIR intervention
const Demography * demography;
const HealthSystem * healthSystem;	// May be replaced by a changeHS intervention
const CaseManagements * caseManagements;	// Optional (may be NULL)
const Parameters * parameters;
//Proteome * proteome;

// Initialized (derived) values:
double parameterValues[Params::MAX];
std::map<int,const Intervention*> timedInterventions;

// Initialization functions:

void initParameterValues() {
  // initialize all to zero
  for (size_t i = 0; i < Params::MAX; ++i)
    parameterValues[i] = 0;
  // set parameters
  const Parameters::ParameterSequence& paramSeq = parameters->getParameter();
  for (Parameters::ParameterConstIterator it = paramSeq.begin(); it != paramSeq.end(); ++it) {
    int i = it->getNumber();
    if (i < 0 || i >= Params::MAX)
      std::cerr << "Warning: parameter with invalid index; ignoring." << std::endl;
    else
      parameterValues[i] = it->getValue();
  }
}

void initTimedInterventions() {
  if (!interventions->getTimed().present()) return;
  const Timed::InterventionSequence& interventionSeq =
      interventions->getTimed().get().getIntervention();
  for (Timed::InterventionConstIterator it (interventionSeq.begin()); it != interventionSeq.end(); ++it) {
    int time = it->getTime();
    if (timedInterventions.count(time)) {
      ostringstream msg;
      msg << "Error: multiple timed interventions with time: " << time;
      throw xml_scenario_error (msg.str());
    }
    timedInterventions[time] = &(*it);
  }
}


void createDocument(std::string lXmlFile) {
  xmlFileName = lXmlFile;
  //Parses the document
    scenario = (parseScenario (lXmlFile)).release();
    if (scenario->getSchemaVersion() < OLDEST_COMPATIBLE) {
      ostringstream msg;
      msg << "Input scenario.xml uses an outdated schema version; please update with SchemaTranslator. Current version: " << SCHEMA_VERSION;
      throw xml_scenario_error (msg.str());
    }
    if (scenario->getSchemaVersion() > SCHEMA_VERSION)
      throw xml_scenario_error ("Error: new schema version unsupported");
    
    monitoring = &scenario->getMonitoring();
    interventions = &scenario->getInterventions();
    entoData = &scenario->getEntoData();
    demography = &scenario->getDemography();
    if (scenario->getHealthSystem().present())
      healthSystem = &scenario->getHealthSystem().get();
    else
      healthSystem = NULL;
    caseManagements = scenario->getCaseManagements().present() ?
        &scenario->getCaseManagements().get() : NULL;
    parameters = &scenario->getParameters();
    //proteome = &scenario->getProteome();
    
    initParameterValues();
    initTimedInterventions();
}

void cleanDocument() {
  if (documentChanged) {
    // get the "basename" (file name without path) of xmlFileName as a C string:
    const char* lastFS = strrchr (xmlFileName.c_str(), '/');
    const char* lastBS = strrchr (xmlFileName.c_str(), '\\');
    const char* baseName = lastBS > lastFS ? lastBS : lastFS;
    if (baseName == NULL)	// no path separator found; use whole string
      baseName = xmlFileName.c_str();
    else
      ++baseName;		// start at next character
    
    ofstream outStream (baseName);
    ostringstream schema;
    schema << "scenario_" << SCHEMA_VERSION << ".xsd";
    
    xml_schema::NamespaceInfomap map;
    map[""].name = "";
    map[""].schema = schema.str();
    serializeScenario (outStream, *scenario, map);
    
    outStream.close();
  }
  
  // Destructors should handle cleanup
  if (scenario != NULL)
    delete scenario;
}

const Monitoring& getMonitoring() {
  return *monitoring;
}
const Interventions& getInterventions() {
  return *interventions;
}
const EntoData& getEntoData() {
  return *entoData;
}
const Demography& getDemography() {
  return *demography;
}
const CaseManagements* getCaseManagements() {
  return caseManagements;
}
const HealthSystem& getHealthSystem() {
  if (healthSystem == NULL)
    throw xml_scenario_error("heathSystem element requested but not present");
  return *healthSystem;
}

void changeHealthSystem (const HealthSystem* hs) {
  healthSystem = hs;
}

double getParameter (size_t i) {
  return parameterValues[i];
}


// ----- Member access functions (bridges) -----
// This is largely unmodified from the old xerces version.

int get_simulation_duration(){ 
  return scenario->getSimulationDuration();
}

double get_detectionlimit(){ 
  return  monitoring->getSurveys().getDetectionLimit();	
}

int get_summary_option(){ 
  return monitoring->getSurveys().getSummaryOption();
}

int get_model_version(){ 
  return scenario->getModelVersion();
}

int get_mode(){
  return scenario->getMode();
}

int get_assim_mode(){ 
  return scenario->getAssimMode();
}

int get_wu_id(){ 
  return scenario->getWuID();
}

double get_maximum_ageyrs(){ 
  return scenario->getMaximumAgeYrs();
}

double get_lowerbound(){ 
  return monitoring->getAgeGroup().getLowerbound();
}

const Intervention* getInterventionByTime(int time){
  std::map<int,const Intervention*>::iterator i = timedInterventions.find (time);
  if (i != timedInterventions.end())
    return i->second;
  else
    return NULL;
}

int get_number_of_surveys(){ 
  return monitoring->getSurveys().getSurveyTime().size();
}

int get_time_of_survey(int index){
  return monitoring->getSurveys().getSurveyTime()[index];
}

int get_analysis_no(){	 
  return scenario->getAnalysisNo();
}

int get_populationsize(){ 
  return scenario->getPopSize();
}


double get_growthrate(){
  if (demography->getGrowthRate().present())
    return demography->getGrowthRate().get();
  else
    return 0.;
}

int get_latentp(){ 
  return parameters->getLatentp();
}

int get_interval(){ 
  return parameters->getInterval();
}

int getISeed(){
  return parameters->getIseed();
}

