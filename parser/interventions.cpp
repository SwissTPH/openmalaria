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

#ifndef INTERVENTIONS_CPP
#define INTERVENTIONS_CPP

#include "Names.h"

#include "vaccineDescription.cpp"
#include "timed.cpp"
#include "ITNDescription.cpp"
#include "IPTDescription.cpp"
#include "continuous.cpp"

#include "MalariaNode.h"
#include <iostream>
#include <math.h>
using namespace std;


class Interventions : public MalariaNode{
private:	
	Timed * timed;
	bool b_timed;

	int numVacineDescription;
	//So far only one IPTI possible
	VaccineDescription** vaccineDescriptions;
	int currentVaccineDescription;
	ITNDescription * iTNDescription;
	bool b_ITNDescription;

	//So far only one IPTI possible
	IPTDescription * iPTDescription;
	bool b_IPTDescription;

	Continuous * continuous;
	bool b_continuous;

public:
	Interventions(DOMNode * node): currentVaccineDescription(0){
		createNode(this,node);
	}

	~Interventions(void) {
		if (b_timed){
			delete timed;
		}
		if (numVacineDescription>0){
			for (int i = 0; i < numVacineDescription; i++){
				delete vaccineDescriptions[i];	}
			delete [] vaccineDescriptions;		
		}

		if (b_ITNDescription){
			delete iTNDescription;}
		if (b_IPTDescription){
			delete iPTDescription;}
		if (b_continuous) {
			delete continuous;}
	}

	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){		
		b_timed = false;
		numVacineDescription = 0;
		for (int i =0, length = list->getLength(); i < length ;i++){
			DOMNode * node = list->item(i);
			if (Converter::equals(node,s_VACCINE_DESCRIPTION)){
				numVacineDescription++;
			}
			vaccineDescriptions = new VaccineDescription*[numVacineDescription];
		}
		b_ITNDescription = false;
		b_IPTDescription = false;
		b_continuous = false;
	}

	void addChild(DOMNode * child){	

		if (Converter::equals(child,s_TIMED)){
			timed =  new Timed(child);			
			b_timed = true;
		}		

		if (Converter::equals(child,s_VACCINE_DESCRIPTION)){
			vaccineDescriptions[currentVaccineDescription++] = new VaccineDescription(child);
		}		

		if (Converter::equals(child,s_ITN_DESCRIPTION)){
			iTNDescription = new ITNDescription(child);
			b_ITNDescription = true;
		}

		if (Converter::equals(child,s_IPT_DESCRIPTION)){
			iPTDescription = new IPTDescription(child);
			b_IPTDescription = true;
		}

		if (Converter::equals(child,s_CONTINUOUS)){
			continuous = new Continuous(child);
			b_continuous = true;
		}	
		//name xml_ipti_description differs from xml_iptiDescription, but the same object is called (IPTDescription)
		if (Converter::equals(child,s_IPTIDESCRIPTION)){
			iPTDescription = new IPTDescription(child);
			b_IPTDescription = true;
		}

	}

	Timed * getTimed(){
		return timed;
	}

	bool isTimed(){
		return b_timed;
	}

	Continuous * getContinuous(){
		return continuous;
	}

	bool isContiuous(){
		return b_continuous;
	}

	bool isVaccineDescription(){
		return numVacineDescription>0;
	}

	int getIndexForType(int type){
		for (int index = 0; index < numVacineDescription; index++){
			if(type==vaccineDescriptions[index]->getVaccineType())
				return index;
		}
		return MISSING_VALUE;
	}

	int getVaccineType(){
		if (numVacineDescription==0){
			return 0;
		}
		int type=0;
		for (int i = 0; i < numVacineDescription; i++){
			type+=(int)pow(2.0,vaccineDescriptions[i]->getVaccineType());	
		}
		return type;
	}

private : VaccineDescription * getVaccineDescriptionByType(int type){
			  return vaccineDescriptions[getIndexForType(type)];
		  }

public : double getEfficacyByType(int type,int dose){
			 if (numVacineDescription==0||getIndexForType(type)==MISSING_VALUE){ 
				 return MISSING_VALUE;}
			 int numberOfDoses = vaccineDescriptions[getIndexForType(type)]->getNumInitialEfficacy();
			 return (dose >= numberOfDoses)?
				 vaccineDescriptions[getIndexForType(type)]->getInitialEfficacy(numberOfDoses-1)->getValue() :
			 vaccineDescriptions[getIndexForType(type)]->getInitialEfficacy(dose)->getValue();
		 }

		 double getHalfLifeByType(int type){
			 if (numVacineDescription==0||getIndexForType(type)==MISSING_VALUE){
				 return MISSING_VALUE;
			 }
			 return  vaccineDescriptions[getIndexForType(type)]->getHalfLifeYears()->getValue();
		 }

		 double getEfficacyBByType(int type){
			 if (numVacineDescription==0||getIndexForType(type)==MISSING_VALUE) {
				 return MISSING_VALUE;
			 }
			 return  vaccineDescriptions[getIndexForType(type)]->getEfficacyB()->getValue();
		 }

		 double getTargetAgeYrs(int index) {
			 if (!b_continuous || continuous->getNumVaccineDoses() == 0) { 
				 return MISSING_VALUE;
			 }	
			 return (index >= continuous->getNumVaccineDoses())?
				 MISSING_VALUE :
			 continuous->getVaccine(index)->getTargetAgeYrs();
		 }

		 double getCoverageEpi(int index) {
			 if (!b_continuous || continuous->getNumVaccineDoses() == 0) { 
				 return MISSING_VALUE;
			 }	
			 return (index >= continuous->getNumVaccineDoses())?
				 MISSING_VALUE :
			 continuous->getVaccine(index)->getCoverage();
		 }

		 int getNumInitEff(){
			 return vaccineDescriptions[0]->getNumInitialEfficacy();
		 }
		
		 int getNumEpiDoses(){
			 return continuous->getNumVaccineDoses();
		 }

		 bool isITNDescription(){
			 return b_ITNDescription;
		 }

		 ITNDescription * getITNDescription(){
			 return iTNDescription;
		 }
		 bool isIPTDescription(){
			 return b_IPTDescription;
		 }

		 IPTDescription * getIPTDescription(){
			 return iPTDescription;
		 }

		 double getIPTiTargetAgeYrs(int index) {
			 if (!b_continuous || continuous->getNumIPTi() == 0) { 
				 return MISSING_VALUE;
			 }	
			 int numberOfIPTiTreatments = continuous->getNumIPTi();
			 return (index >= numberOfIPTiTreatments)?
				 continuous->getIPTi(numberOfIPTiTreatments-1)->getTargetAgeYrs() :
			 continuous->getIPTi(index)->getTargetAgeYrs();
		 }


#ifdef _LOG
		 void debug(){}
#endif

};

#endif
