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

// Other parameter-getters
// WARNING: not all these functions may exist

double get_parameter(int); 

double get_detectionlimit(); 
double get_eir_daily(int); 
Anopheles* getAnopheles(std::string name);
int is_survey(int); 
int get_number_of_surveys(); 
int get_time_of_survey(int ); 
int get_summary_option();	 
int get_model_version(); 
int get_mode(); 
int get_assim_mode(); 
int get_wu_id();
int get_release();
double get_maximum_ageyrs(); 
double get_lowerbound(); 
int get_number_of_agegroups(); 
double get_upperbound(int ); 
int get_analysis_no(); 
int get_populationsize(); 
double get_growthrate(); 

int get_simulation_duration(); 
double get_p_sequelae(int ); 
int get_number_of_cfrgroups(); 
double get_cfr_lb(int ); 
double get_cfr(int ); 

int get_vaccine_type(); 
int get_vaccine_effect(); 
int get_number_of_epi_doses(); 
int get_number_of_init_eff(); 
double get_target_age_yrs(int ); 
double get_efficacy(int ,int ); 
double get_efficacy_b(int ); 
double get_vaccine_halflife_yrs(int ); 
double get_coverage_epi_vaccine(int ); 
double get_pu0(); 
double get_pu1(); 
double get_sporogony_gonotrophy(); 
double get_itn_halflife_yrs(); 
double get_demo_lowerbound(); 
double get_demo_upperbound(int ); 
double get_popperc(int); 

double get_delta(); 
int get_latentp(); 
int get_interval(); 
double get_ipti_effect(); 
double get_genotype_freq(int ); 
double get_genotype_acr(int ); 
int get_genotype_proph(int ); 
int get_genotype_tolperiod(int ); 
double get_genotype_atten(int ); 
int get_is_ipti(); 
int get_number_of_genotypes(); 
double get_ipti_coverage(int ); 
double get_ipti_target_age_yrs(int ); 
int get_number_of_ipti_doses(); 

int get_number_of_proteins(); 
char* get_protein_name_i(int); 
void get_protein_name(char*, int, int); 
int get_protein_number_of_mutations(int); 
int get_protein_mutation_position(int, int); 
int get_number_of_proteome_instances(); 
double get_pi_proportion(int); 
double get_pi_fitness(int); 
int get_pi_number_of_alleles(int); 
char* get_allele_name_i(int, int); 
void get_allele_name(char*, int, int, int); 
int get_allele_cnv(int, int); 
char* get_allele_aminos_i(int, int); 
void get_allele_aminos(char*, int, int, int); 

int get_decision_id(int, double); 
int get_n_medicate(int, double); 
double get_cmp_qty(int, int, double); 
int get_cmp_time(int, int, double); 
void get_cmp_name(char*, int, int, int, double); 
int getISeed(); 

#endif

