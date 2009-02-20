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

/**
* C++ wrapper for the xml scenario
*/
#ifndef SCENARIO_CPP
#define SCENARIO_CPP

#include "Names.h"

#include "monitoring.cpp"
#include "converter.cpp"
#include "demography.cpp"
#include "interventions.cpp"
#include "entoData.cpp"
#include "healthSystem.cpp"
#include "caseManagements.cpp"
#include "parameters.cpp"
#include "proteome.cpp"

#include <xercesc/dom/DOM.hpp>
#include "MalariaNode.h"
XERCES_CPP_NAMESPACE_USE

#include <iostream>

using namespace std;



class Scenario : public MalariaNode
{
private:
	int modelVersion;
	int mode;
	int analysisNo;
	char * name;
	double maximumAgeYears;
	int popSize;
	int simulationDuration;
	int wuID;
	int assimMode;
	double growthrate;
  int release;

	//Some booleans to know if the values have been defined
	bool b_demography;
	bool b_monitoring;
	bool b_interventions;
	bool b_entoData;
	bool b_healthSystem;
	bool b_caseManagements;
	bool b_parameters;
	bool b_modelVersion;
	bool b_mode;
	bool b_wu_id;
	bool b_assim_mode;
  bool b_release;
  bool b_growthrate;
	bool b_proteome;

	Demography * demography;
	Monitoring * monitoring;
	Interventions * interventions;
	EntoData *entoData;
	HealthSystem * healthSystem;
	CaseManagements * caseManagements;
	Parameters * parameters;
    Proteome * proteome;
	
	
public:
	Scenario(DOMNode * node){			
		proteome = NULL;
		createNode(this,node);
	}
	
	~Scenario(void)	{
		delete [] name;
		//We delete the object only if they were defined
		if (b_demography) delete demography;
		if (b_monitoring) delete monitoring;
		if (b_interventions) delete interventions;
		if (b_entoData) delete entoData;
		if (b_healthSystem) delete healthSystem;
		if (b_caseManagements) delete caseManagements;
		if (b_parameters) delete parameters;
		if (b_proteome) delete proteome;
	}


	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		//Initialise the booleans
		b_demography = false;
		b_monitoring = false;
		b_interventions = false;
		b_entoData = false;
		b_healthSystem = false;
		b_caseManagements = false;
		b_parameters = false;
		b_modelVersion = false;
		b_mode = false;
		b_wu_id = false;
		b_assim_mode = false;
    b_release = false;
		b_growthrate = false;
		b_proteome = false;

		//Get the different values
		maximumAgeYears = Converter::parseDouble(s_MAXIMUM_AGE_YEARS,map);
		analysisNo = Converter::parseInt(s_ANALYSIS_NO,map);
		name = Converter::getValue(s_NAME,map);	
		popSize = Converter::parseInt(s_POP_SIZE,map);
		simulationDuration = Converter::parseInt(s_SIMULATION_DURATION,map);
		if (Converter::contains(s_MODEL_VERSION,map)){
			modelVersion = Converter::parseInt(s_MODEL_VERSION,map);
			b_modelVersion=true;
		}
		if (Converter::contains(s_MODE,map)){
			mode = Converter::parseInt(s_MODE,map);
			b_mode=true;
		}
		if (Converter::contains(s_ASSIM_MODE,map)){
			assimMode = Converter::parseInt(s_ASSIM_MODE,map);
			b_assim_mode=true;
		}
    if (Converter::contains(s_RELEASE,map)){
      release = Converter::parseInt(s_RELEASE,map);
			b_release=true;
		}
		if (Converter::contains(s_WU_ID,map)){
			wuID = Converter::parseInt(s_WU_ID,map);
			b_wu_id=true;
		}
		if (Converter::contains(s_GROWTHRATE,map)){
			growthrate = Converter::parseDouble(s_GROWTHRATE,map);
			b_growthrate=true;
		}

	}

	void addChild(DOMNode * child){
		//Get the demography
		if (Converter::equals(child,s_DEMOGRAPHY)){
			b_demography = true;
			demography =  new Demography(child);		
			return;
		}

		//Get the monitoring
		if (Converter::equals(child,s_MONITORING)){
			b_monitoring = true;
			monitoring =  new Monitoring(child);			
			return;
		}
		//Get the interventions
		if (Converter::equals(child,s_INTERVENTIONS)){
			b_interventions = true;
			interventions =  new Interventions(child);			
			return;
		}

		//Get the entomolgy
		if (Converter::equals(child,s_ENTO_DATA)){
			b_entoData = true;
			entoData =  new EntoData(child);			
			return;
		}
		//Get the healthSystem
		if (Converter::equals(child,s_HEALTH_SYSTEM)){
			b_healthSystem = true;
			healthSystem =  new HealthSystem(child);			
			return;
		}
		//Get the caseManagements 
		if (Converter::equals(child,s_CASE_MANAGEMENTS)){
			b_caseManagements = true;
			caseManagements = new CaseManagements(child);			
			return;
		}
		//Get the parameters
		if (Converter::equals(child,s_PARAMETERS)){
			b_parameters = true;
			parameters = new Parameters(child);
			return;
		}
		//Get the proteome
		if (Converter::equals(child,s_PROTEOME)){
			b_proteome = true;
			proteome = new Proteome(child);
			return;
		}
	}

	Monitoring * getMonitoring(){
		return monitoring;
	}

	EntoData * getEntomology(){
		return entoData;
	}

	Demography * getDemography(){
		return demography;
	}

	HealthSystem * getHealthSystem(){
		return healthSystem;
	}

	CaseManagements * getCaseManagements(){
		return caseManagements;
	}
	//
	// Return the parameters
	//
	Parameters * getParameters(){
		return parameters;
	}

	Proteome * getProteome(){
		return proteome;
	}

	int getSimulationDuration(){
		return simulationDuration;
	}

	int getMode(){
		if(b_mode){
		return mode;
		}
		else{
			return MISSING_VALUE;
		}
	}

	int getAssimMode(){
		if(b_assim_mode){
		return assimMode;
		}
		else{
			return MISSING_VALUE;
		}
	}

  int getRelease(){
		if(b_release){
      return release;
		}
		else{
			return MISSING_VALUE;
		}
	}

	int getWuID(){
		if(b_wu_id){
		return wuID;
		}
		else{
			return MISSING_VALUE;
		}
	}
	
	int getModelVersion(){
		if(b_modelVersion){
		return modelVersion;
		}
		else{
			return MISSING_VALUE;
		}
	}

	double getMaximumAgeYears(){
		return maximumAgeYears;
	}

	int getAnalysisNo(){
		return analysisNo;
	}

	Interventions * getInterventions(){
		return interventions;
	}

	int getPopulationSize(){
		return popSize;
	}

	double getGrowthrate(){
		if(b_growthrate){
			return growthrate;
		}
		else{
			//return zero as default
			return 0.0;
		}
	}


#ifdef _LOG
	void Scenario::debug(){
		cerr << "<Scenario" << "\t" <<
			s_ANALYSIS_NO << " " << analysisNo << "\t" <<
			s_NAME << " " << name << "\t" <<
			s_POP_SIZE << " " << popSize << "\t" <<
			s_SIMULATION_DURATION << " " << simulationDuration << "\t" <<
			">\n";	
	}
#endif
};

#endif
