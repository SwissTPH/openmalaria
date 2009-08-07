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
// Naming convention can be changed, but was chosen to minimize changes in the
// code below.

int get_simulation_duration(){ 
  return scenario->getSimulationDuration();
}

double get_detectionlimit(){ 
  return  monitoring->getSurveys().getDetectionLimit();	
}

int is_survey(int time) {
  const Surveys::SurveyTimeSequence& times = monitoring->getSurveys().getSurveyTime();
  for (Surveys::SurveyTimeConstIterator it = times.begin(); it != times.end(); ++it)
    if (time == *it)
      return true;
  return false;
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


double get_demo_lowerbound(){ 
  return demography->getAgeGroup().getLowerbound();
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
double get_delta(){ 
  return parameters->getDelta();
}

/*
int get_number_of_proteins(){ 
  return proteome ? proteome->getContent().getNumProteins() : 0;
}

void get_fortran_string(char* result, int length, char* orig) {
  int gotNull = 0;
  for(int i=0;i<length; i++) {
    if (orig[i] == 0) {
      gotNull = 1;
    }
    if (gotNull) {
      result[i] = ' ';
    }
    else {
      result[i] = orig[i];
    }
  }
}

char* get_protein_name_i(int index){
  return proteome->getContent().getProtein(index - 1).getName();
}

void get_protein_name(char* result, int length, int index) { 
  get_fortran_string(
                     result, length,
                     proteome->getContent().getProtein(index - 1).getName());
}

int get_protein_number_of_mutations(int index){ 
  return proteome->getContent().getProtein(index - 1).getNumMutations();
}

int get_protein_mutation_position(int gindex, int mindex){ 
  return proteome->getContent()->getProtein(gindex - 1)->getMutation(mindex - 1)->getPosition();
}

int get_number_of_proteome_instances(){
  return proteome ? proteome->getDistribution()->getNumProteomeInstances() : 0;
}

double get_pi_proportion(int giindex){ 
  return proteome->getDistribution()->getProteomeInstance(giindex - 1)->getProportion();
}

double get_pi_fitness(int giindex){ 
  return proteome->getDistribution()->getProteomeInstance(giindex - 1)->getFitness();
}

/* Probably to deprecate *//*
int get_pi_number_of_alleles(int giindex){ 
  return proteome->getDistribution()->getProteomeInstance(giindex - 1)->getNumAlleles();
}

char* get_allele_name_i(int giindex, int aindex){
  return proteome->getDistribution()->getProteomeInstance(giindex - 1)->getAllele(aindex - 1)->getName();
}

void get_allele_name(char* result, int length, int giindex, int aindex) { 
  get_fortran_string(
                     result, length,
                     proteome->getDistribution()->getProteomeInstance(giindex - 1)->getAllele(aindex - 1)->getName());
}

/* For future use *//*
int get_allele_cnv(int giindex, int aindex){ 
  return proteome->getDistribution()->getProteomeInstance(giindex - 1)->getAllele(aindex - 1)->getCNV();
}

char* get_allele_aminos_i(int giindex, int aindex) {
  return proteome->getDistribution()->getProteomeInstance(giindex - 1)->getAllele(aindex - 1)->getAminos();
}

void get_allele_aminos(char* result, int length, int giindex, int aindex) { 
  get_fortran_string(
                     result, length,
                     proteome->getDistribution()->getProteomeInstance(giindex - 1)->getAllele(aindex - 1)->getAminos());
}*/

/*
double get_cmp_qty(int decisionID, int medicateID, double age){ 
  return caseManagements->getCaseManagementbyAge(age)->getDecisions()->getDecision(decisionID-1)->getMedicate(medicateID-1)->getQTY();
}

int get_cmp_time(int decisionID, int medicateID, double age){ 
  return caseManagements->getCaseManagementbyAge(age)->getDecisions()->getDecision(decisionID-1)->getMedicate(medicateID-1)->getTime();
}

void get_cmp_name(char* result, int length,int decisionID, int medicateID, double age){ 
  get_fortran_string(
                     result, length,
                     caseManagements->getCaseManagementbyAge(age)->getDecisions()->getDecision(decisionID-1)->getMedicate(medicateID-1)->getName());
}*/

int getISeed(){
  return parameters->getIseed();
}

