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

using namespace std;
#include "inputData.h"
#include <math.h>
#include <iostream>

#include "document.cpp"
#include "anopheles.cpp"
#include "DOMParser.h"

Document * document;
Scenario * scenario ;
Monitoring * monitoring;
Interventions * interventions;
EntoData * entoData;
Demography * demography;
HealthSystem * healthSystem;
CaseManagements * caseManagements;
Parameters * parameters;
Proteome * proteome;

double cureRate [3];
double parasitesCleared [3];

void initialiseCureRate(){
	//getThe first line of the drug regimen
	char *firstLineDrug = healthSystem->getDrugRegimen()->getFirstLine();
	//We get the ACR depending on the name of firstLineDrug.
	double curerateFirstLine = healthSystem->getInitialACR()->getACRByName(firstLineDrug);
	double pSeekOfficialCareUncomplicated1 = healthSystem->getPSeekOfficialCareUncomplicated1()->getValue();
	double cureRateSelfTreatment = healthSystem->getInitialACR()->getSelfTreatment()->getValue();
	double pSelfTreatment = healthSystem->getPSelfTreatUncomplicated()->getValue();
	//Calculate curerate 0
    if ((pSeekOfficialCareUncomplicated1+pSelfTreatment)>0){
        cureRate[0] = (curerateFirstLine*pSeekOfficialCareUncomplicated1+
        cureRateSelfTreatment*pSelfTreatment)/
               (pSeekOfficialCareUncomplicated1+pSelfTreatment);
    }else{
         cureRate[0]=curerateFirstLine;
    }

	//Calculate curerate 1
	char * secondLineDrug = healthSystem->getDrugRegimen()->getSecondLine();	
	double cureRateSecondLine = healthSystem->getInitialACR()->getACRByName(secondLineDrug);
	cureRate[1] = cureRateSecondLine;

	//Calculate curerate 2
	char * inpatient = healthSystem->getDrugRegimen()->getInpatient();
	double cureRateInpatient = healthSystem->getInitialACR()->getACRByName(inpatient);
	cureRate[2] = cureRateInpatient;	
}

void initialiseParasitesCleared(){
	char *firstLineDrug = healthSystem->getDrugRegimen()->getFirstLine();
	char * secondLineDrug = healthSystem->getDrugRegimen()->getSecondLine();	
	double pSeekOfficialCareUncomplicated1 = healthSystem->getPSeekOfficialCareUncomplicated1()->getValue();

	double complianceFirstLine = healthSystem->getCompliance()->getACRByName(firstLineDrug);
	double complianceSecondLine = healthSystem->getCompliance()->getACRByName(secondLineDrug);

	double cureRateFirstLine = healthSystem->getInitialACR()->getACRByName(firstLineDrug);
	double cureRateSecondLine = healthSystem->getInitialACR()->getACRByName(secondLineDrug);

	double nonCompliersEffectiveFirstLine = healthSystem->getNonCompliersEffective()->getACRByName(firstLineDrug);
	double nonCompliersEffectiveSecondLine = healthSystem->getNonCompliersEffective()->getACRByName(secondLineDrug);

	double pSelfTreatment = healthSystem->getPSelfTreatUncomplicated()->getValue();
	double complianceSelfTreatment = healthSystem->getCompliance()->getSelfTreatment()->getValue();
	double cureRateSelfTreatment = healthSystem->getInitialACR()->getSelfTreatment()->getValue();
	//calculate parasitesCleared 0
	if ((pSeekOfficialCareUncomplicated1+pSelfTreatment)>0){
	parasitesCleared[0] = (pSeekOfficialCareUncomplicated1
		* (complianceFirstLine * cureRateFirstLine + (1 - complianceFirstLine)
		* nonCompliersEffectiveFirstLine) + pSelfTreatment
		* (complianceSelfTreatment * cureRateSelfTreatment + (1 -
		complianceSelfTreatment)
		*nonCompliersEffectiveFirstLine))
		/( pSeekOfficialCareUncomplicated1 + pSelfTreatment);
	}else{
	parasitesCleared[0] = 0;
	}
	//calculate parasitesCleared 1
	parasitesCleared[1] = complianceSecondLine * cureRateSecondLine
		+ (1-complianceSecondLine)
		* nonCompliersEffectiveSecondLine;

	//calculate parasitesCleared 2 : cool :)
	parasitesCleared[2] = 0;
}


/**
*	Reads the document in the xmlFile
* Return if the result was a success
*/
bool createDocument(const char * lXmlFile){
	//Parses the document
	try{
		document = (Document*) parseMalaria(lXmlFile);
		//Get the different root nodes.
		scenario = document->getScenario();
		monitoring = scenario->getMonitoring();
		entoData = scenario->getEntomology();
		interventions = scenario->getInterventions();
		demography = scenario->getDemography();
		caseManagements = scenario->getCaseManagements();
		healthSystem = scenario->getHealthSystem();
		parameters = scenario->getParameters();
                proteome = scenario->getProteome();

		//Initialises the different values of the health system
		initialiseCureRate();
		initialiseParasitesCleared();
	}
	catch (...){
		cerr << "IDB\n";
		return false;
	}
	return true;
}

void cleanDocument(){
	delete document;
	//all intern classe are deleted too.
}

int get_simulation_duration(){ 
  return scenario->getSimulationDuration();
}

double get_detectionlimit(){ 
  return  monitoring->getSurveys()->getDetectionLimit();	
}

int is_survey(int time){ 
  return monitoring->isSurvey(time);
}

int get_summary_option(){ 
  return monitoring->getSurveys()->getSummaryOption();
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

int get_release(){ 
  return scenario->getRelease();
}

int get_wu_id(){ 
  return scenario->getWuID();
}

double get_maximum_ageyrs(){ 
  return scenario->getMaximumAgeYears();
}

double get_lowerbound(){ 
  return monitoring->getAgeGroup()->getLowerBound();
}

int get_number_of_agegroups(){ 
  return monitoring->getAgeGroup()->getNumGroups();
}

double get_upperbound(int index){ 
  return monitoring->getAgeGroup()->getGroup(index)->getUpperBound();	
}

int get_intervention(int time){			 
	if (!interventions->isTimed()) 
		return NO_INTERVENTION;
	Intervention* intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL)
		return NO_INTERVENTION;
	int interventionCode = NO_INTERVENTION;
	if (intervention->isIRS())
		interventionCode += (int) pow(2.0,IRS_INTERVENTION);
	if (intervention->isMDA())
		interventionCode += (int) pow(2.0,MDA_INTERVENTION);
	if (intervention->isIPTi())
		interventionCode += (int) pow(2.0,IPTi_INTERVENTION);
	if (intervention->isVaccine())
		interventionCode += (int) pow(2.0,VACCINE_INTERVENTION);
	if (intervention->isChangeEIR()){
		interventionCode += (int) pow(2.0,CHANGE_EIR_INTERVENTION);
		delete entoData;
		entoData=intervention->getChangeEIR();
	}
	if (intervention->isChangeHS()){
		interventionCode += (int) pow(2.0,CHANGE_HS_INTERVENTION);
		delete healthSystem;
		healthSystem=intervention->getChangeHS();
		//Initialises the different values of the health system
		initialiseCureRate();
		initialiseParasitesCleared();
	}
	return interventionCode;
}

double get_maxage_mda(int time){ 
	if (!interventions->isTimed()) 
		return MISSING_VALUE;
	Intervention * intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL || !intervention->isMDA())
		return MISSING_VALUE;
	Mass * mda = intervention->getMDA();
	return mda->getMaxAge();	
}

double get_minage_mda(int time){ 
	if (!interventions->isTimed()) 
		return MISSING_VALUE;
	Intervention * intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL || !intervention->isMDA())
		return MISSING_VALUE;
	Mass * mda = intervention->getMDA();
	return mda->getMinAge();	
}

double get_coverage_mda(int time){ 
	if (!interventions->isTimed()) 
		return MISSING_VALUE;
	Intervention * intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL || !intervention->isMDA())
		return MISSING_VALUE;
	Mass * mda = intervention->getMDA();
	return mda->getCoverage();	
}

double get_maxage_vaccine(int time){ 
	if (!interventions->isTimed()) 
		return MISSING_VALUE;
	Intervention * intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL || !intervention->isVaccine())
		return MISSING_VALUE;
	Mass * vaccine = intervention->getVaccine();
	return vaccine->getMaxAge();	
}

double get_minage_vaccine(int time){ 
	if (!interventions->isTimed()) 
		return MISSING_VALUE;
	Intervention * intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL || !intervention->isVaccine())
		return MISSING_VALUE;
	Mass * vaccine = intervention->getVaccine();
	return vaccine->getMinAge();	
}

double get_maxage_ipti(int time){ 
	if (!interventions->isTimed()) 
		return MISSING_VALUE;
	Intervention * intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL || !intervention->isIPTi())
		return MISSING_VALUE;
	Mass * ipti = intervention->getIPTi();
	return ipti->getMaxAge();	
}

double get_minage_ipti(int time){ 
	if (!interventions->isTimed()) 
		return MISSING_VALUE;
	Intervention * intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL || !intervention->isIPTi())
		return MISSING_VALUE;
	Mass * ipti = intervention->getIPTi();
	return ipti->getMinAge();	
}

double get_coverage_ipti(int time){ 
	if (!interventions->isTimed()) 
		return MISSING_VALUE;
	Intervention * intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL || !intervention->isIPTi())
		return MISSING_VALUE;
	Mass * ipti = intervention->getIPTi();
	return ipti->getCoverage();	
}

double get_coverage_epi_vaccine(int index){ 
	return interventions->getCoverageEpi(index);
}

double get_coverage_mass_vaccine(int time){ 
	if (!interventions->isTimed()) 
		return MISSING_VALUE;
	Intervention * intervention = interventions->getTimed()->getInterventionByTime(time);
	if (intervention == NULL || !intervention->isVaccine())
		return MISSING_VALUE;
	Mass * vaccine = intervention->getVaccine();
	return vaccine->getCoverage();
}

int get_number_of_surveys(){ 
	return monitoring->getSurveys()->getNumSurveys();
}

int get_time_of_survey(int index){ 
	return monitoring->getSurveys()->getSurvey(index);
}

int get_analysis_no(){	 
	return scenario->getAnalysisNo();
}

int get_populationsize(){ 
	return scenario->getPopulationSize();
}


double get_p_gets_treatment(int index){ 

  switch (index){
  case 0:	//treatments[0]
    return healthSystem->getPSeekOfficialCareUncomplicated1()->getValue() + 
    healthSystem->getPSelfTreatUncomplicated()->getValue();

  case 1:	//treatments[1]
    return healthSystem->getPSeekOfficialCareUncomplicated2()->getValue();

  case 2:	//treatments[2]
    return healthSystem->getPSeekOfficialCareSevere()->getValue();

  default :
    cerr <<"No such treatment\n";
    exit(-1);
  }
  
}

double get_curerate(int index){	 
  return cureRate[index];
}

double get_p_parasites_cleared(int index){ 
  return parasitesCleared[index];
}

double get_p_sequelae(int index){		 
    
  //TODO: very simple way to model sequelae, based on only 2 agegroups
  if (0==index){
    return healthSystem->getPSequelaeInpatient()->getByAge(1);
  }
  else{
    return healthSystem->getPSequelaeInpatient()->getByAge(10);
  }
}

int get_number_of_cfrgroups(){ 
  return healthSystem->getCFR()->getNumGroups();
}

double get_cfr(int index){	 
  return healthSystem->getCFR()->getGroup(index)->getCFR();
}

double get_cfr_lb(int index){	 
  return healthSystem->getCFR()->getGroup(index)->getLowerBound();
}

int get_health_system_memory(){ 
  return healthSystem->getHealthSystemMemory();
}

int get_vaccine_type(){ 
  return (interventions->isVaccineDescription())? interventions->getVaccineType() : 0;
}

int get_number_of_epi_doses(){ 
  return interventions->getNumEpiDoses();
}

int get_number_of_init_eff(){ 
  return interventions->getNumInitEff();
}

double get_vaccine_halflife_yrs(int type){ 
  return (interventions->isVaccineDescription())?interventions->getHalfLifeByType(type): MISSING_VALUE;
}

double get_efficacy_b(int type){ 
  return (interventions->isVaccineDescription())?interventions->getEfficacyBByType(type): MISSING_VALUE;
}

double get_target_age_yrs(int index){ 
  return interventions->getTargetAgeYrs(index);
}

double get_efficacy(int type,int index){ 
  return interventions->getEfficacyByType(type,index);
}

double get_pu0(){ 
  return (!interventions->isITNDescription())? 0 : interventions->getITNDescription()->getPu0()->getValue();
}

double get_pu1(){ 
  return (!interventions->isITNDescription())? 0 : interventions->getITNDescription()->getPu1()->getValue();
}

double get_sporogony_gonotrophy(){	 
  return (!interventions->isITNDescription())? 0 : 
		interventions->getITNDescription()->getSporogonyGonotrophy()->getValue();
}

double get_itn_halflife_yrs(){	 
  return (!interventions->isITNDescription())? 
		0 : interventions->getITNDescription()->getHalfLifeYrs()->getValue();
}

double get_demo_lowerbound(){ 
  return demography->getAgeGroup()->getLowerBound();
}

double get_demo_upperbound(int index){ 
  return demography->getAgeGroup()->getGroup(index)->getUpperBound();	
}

double get_popperc(int index){ 
  return demography->getAgeGroup()->getGroup(index)->getPopPercent();	
}

double get_growthrate(){ 
  return scenario->getGrowthrate();
}

double get_eir_daily(int time){ 
  return entoData->getEIRDaily(time);
}

Anopheles * getAnopheles(string name){
  return entoData->getAnopheles(name);
}

double get_parameter(int index){ 
  return parameters->getParameter(index-1)->getValue();	
}

int get_num_parameters(){ 
  return parameters->getNumParameters();
}

double get_latentp(){ 
  return parameters->getLatentp();
}

double get_nspore(){ 
  return parameters->getNspore();
}

int get_interval(){ 
  return parameters->getInterval();
}
double get_delta(){ 
  return parameters->getDelta();
}

double get_ipti_effect(){ 
  return interventions->getIPTDescription()->getIptiEffect();
}

double get_genotype_freq(int index){ 
  return interventions->getIPTDescription()->getGenoType(index-1)->getFreq();
}

double get_genotype_acr(int index){ 
  return interventions->getIPTDescription()->getGenoType(index-1)->getACR();
}

int get_genotype_proph(int index){ 
  return interventions->getIPTDescription()->getGenoType(index-1)->getProph();
}

int get_genotype_tolperiod(int index){ 
  return interventions->getIPTDescription()->getGenoType(index-1)->getTolperiod();
}
double get_genotype_atten(int index){ 
  return interventions->getIPTDescription()->getGenoType(index-1)->getAtten();
}

int get_is_ipti(){ 
  if(interventions->isIPTDescription()){
    return 1;
  }
  else{
    return 0;
  }
}

int get_number_of_genotypes(){ 
  return interventions->getIPTDescription()->getnumGenoTypes();
}

double get_ipti_coverage(int index){ 
  //only call this if you are sure that ->isIPTDescription() = true
  if (!interventions->isContiuous() || !interventions->isIPTDescription()) {
    return MISSING_VALUE;
  }
  int numberOfIPTiTreatments = interventions->getContinuous()->getNumIPTi();
    return (index >= numberOfIPTiTreatments)?
		interventions->getContinuous()->getIPTi(numberOfIPTiTreatments)->getCoverage() :
	interventions->getContinuous()->getIPTi(index)->getCoverage();
}

double get_ipti_target_age_yrs(int index){ 
  return interventions->getIPTiTargetAgeYrs(index);
}

int get_number_of_ipti_doses(){ 
  if (interventions->isContiuous() && interventions->isIPTDescription()){
	return interventions->getContinuous()->getNumIPTi();
  }
  else{
    return 0;
  }	
}

int get_number_of_proteins(){ 
  return proteome ? proteome->getContent()->getNumProteins() : 0;
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
  return proteome->getContent()->getProtein(index - 1)->getName();
}

void get_protein_name(char* result, int length, int index) { 
    get_fortran_string(
		result, length,
		proteome->getContent()->getProtein(index - 1)->getName());
}

int get_protein_number_of_mutations(int index){ 
  return proteome->getContent()->getProtein(index - 1)->getNumMutations();
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

/* Probably to deprecate */
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

/* For future use */
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
}

int get_decision_id(int entryPointID, double age){ 
  switch (entryPointID) {

  case 1 : 
    return caseManagements->getCaseManagementbyAge(age)->getUncomplicatedFirst()->getDecisionID();
    break;
  case 2 : 
    return caseManagements->getCaseManagementbyAge(age)->getUncomplicatedSecond()->getDecisionID();
    break;
  case 3 :
    return caseManagements->getCaseManagementbyAge(age)->getSevere()->getDecisionID();
    break;
  case 4 :
    return caseManagements->getCaseManagementbyAge(age)->getNMF()->getDecisionID();
    break;

  default :
    cerr <<"No such entrypoint\n";
    exit(-1);
  }
}

int get_n_medicate(int decisionID, double age){
  return caseManagements->getCaseManagementbyAge(age)->getDecisions()->getDecision(decisionID-1)->getNumMedicates();
}

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
}

int getISeed(){
  return parameters->getISeed();
}

