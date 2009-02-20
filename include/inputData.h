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
#include <string>
class Anopheles;

/**
* return if the function was a success. The goal of this function is to parse. 
Do not wait to have a policy saying you were the error is. In the best case, the function
will crash and you will be able to use the debugger to know where the error is.
*/
bool createDocument(const char *);

/**
* Some elements in memory have been created. This function deletes the object in memory
*/
void cleanDocument();

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
int get_intervention(int ); 
double get_maxage_mda(int ); 
double get_minage_mda(int ); 
double get_coverage_mda(int ); 
double get_maxage_ipti(int ); 
double get_minage_ipti(int ); 
double get_coverage_ipti(int ); 
double get_maxage_vaccine(int ); 
double get_minage_vaccine(int ); 
double get_coverage_mass_vaccine(int ); 
int get_analysis_no(); 
int get_populationsize(); 
double get_growthrate(); 

int get_simulation_duration(); 
double get_p_gets_treatment(int ); 
double get_curerate(int ); 
double get_p_parasites_cleared(int ); 
double get_p_sequelae(int ); 
int get_number_of_cfrgroups(); 
double get_cfr_lb(int ); 
double get_cfr(int ); 

int get_health_system_memory(); 
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

double get_parameter(int); 
int get_num_parameters(); 
double get_delta(); 
double get_latentp(); 
double get_nspore(); 
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

