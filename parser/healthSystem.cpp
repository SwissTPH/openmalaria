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



#ifndef HEALTH_SYSTEM_CPP
#define HEALTH_SYSTEM_CPP

#include "converter.cpp"
#include "Names.h"
#include "drugRegimen.cpp"
#include "treatmentDetails.cpp"
#include "value.cpp"
#include "byAgeItems.cpp"
#include "CFR.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class HealthSystem : public MalariaNode
{
private:	
	char * name;
	int healthSystemMemory;	
	DrugRegimen * drugRegimen;
	TreatmentDetails * initialACR;
	TreatmentDetails * compliance;
	TreatmentDetails * nonCompliersEffective;
	Value * pSeekOfficialCareUncomplicated1;
	Value * pSeekOfficialCareUncomplicated2;
	Value * pSelfTreatUncomplicated;
	Value * pSeekOfficialCareSevere;
	ByAgeItems * pSequelaeInpatient;
	CFR * cfr;
public:
	HealthSystem(DOMNode * node){
		createNode(this,node);
	};

	~HealthSystem(void) {
		delete name;
		delete drugRegimen;
		delete initialACR;
		delete compliance;
		delete nonCompliersEffective;
		delete pSeekOfficialCareUncomplicated1;
		delete pSelfTreatUncomplicated;
		delete pSeekOfficialCareUncomplicated2;
		delete pSeekOfficialCareSevere;
		delete cfr;
		delete pSequelaeInpatient;
	};
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		name = Converter::getValue(s_NAME,map);
		healthSystemMemory = Converter::parseInt(s_HEALTH_SYSTEM_MEMORY,map);
	}

	void addChild(DOMNode * child){	
		if (Converter::equals(child,s_DRUG_REGIMEN)){
			drugRegimen =  new DrugRegimen(child);			
		}

		if (Converter::equals(child,s_INITIAL_ACR)){
			initialACR = new TreatmentDetails(s_INITIAL_ACR,child);
		}

		if (Converter::equals(child,s_COMPLIANCE)){
			compliance = new TreatmentDetails(s_COMPLIANCE,child);
		}

		if (Converter::equals(child,s_NON_COMPLIERS_EFFECTIVE)){
			nonCompliersEffective = new TreatmentDetails(s_NON_COMPLIERS_EFFECTIVE,child);
		}

		if (Converter::equals(child,s_P_SEEK_OFFICIAL_CARE_UNCOMPLICATED_1)){
			pSeekOfficialCareUncomplicated1 = 
				new Value(s_P_SEEK_OFFICIAL_CARE_UNCOMPLICATED_1,child);
		}

		if (Converter::equals(child,s_P_SEEK_OFFICIAL_CARE_UNCOMPLICATED_2)){
			pSeekOfficialCareUncomplicated2 = 
				new Value(s_P_SEEK_OFFICIAL_CARE_UNCOMPLICATED_2,child);
		}

		if (Converter::equals(child,s_P_SELF_TREAT_UNCOMPLICATED)){
			pSelfTreatUncomplicated = 
				new Value(s_P_SELF_TREAT_UNCOMPLICATED,child);
		}

		if (Converter::equals(child,s_P_SEEK_OFFICIAL_CARE_SEVERE)){
			pSeekOfficialCareSevere = 
				new Value(s_P_SEEK_OFFICIAL_CARE_SEVERE,child);
		}
		
		if (Converter::equals(child,s_P_SEQUELAE_INPATIENT)){
			pSequelaeInpatient = new ByAgeItems(child);
		}

		if (Converter::equals(child,s_CFR)){
			cfr = new CFR(child);
		}
	}

	Value * getPSeekOfficialCareUncomplicated1(){
		return pSeekOfficialCareUncomplicated1;
	}

	Value * getPSeekOfficialCareUncomplicated2(){
		return pSeekOfficialCareUncomplicated2;
	}

	Value * getPSelfTreatUncomplicated(){
		return pSelfTreatUncomplicated;
	}

	Value * getPSeekOfficialCareSevere(){
		return pSeekOfficialCareSevere;
	}

	TreatmentDetails * getInitialACR(){
		return initialACR;
	}

	TreatmentDetails * getCompliance(){
		return compliance;
	}
	
	TreatmentDetails * getNonCompliersEffective(){
		return nonCompliersEffective;
	}

	DrugRegimen * getDrugRegimen(){
		return drugRegimen;
	}	

	int getHealthSystemMemory(){
		return healthSystemMemory;
	}

	CFR * getCFR(){
		return cfr;
	}

	ByAgeItems * getPSequelaeInpatient(){
		return pSequelaeInpatient;
	}

#ifdef _LOG
	void debug(){
		cerr << "<healthSystem " <<
			"\tname " << name <<
			"\thealthSystemMemory " << healthSystemMemory <<
			"\t>\n";
	}
#endif

};

#endif
