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

#ifndef CASE_MANAGEMENT_CPP
#define CASE_MANAGEMENT_CPP

#include "Names.h"
#include "converter.cpp"
#include "entryPoint.cpp"
#include "decisions.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class CaseManagement : public MalariaNode{
private:	
	//Uncomplicated episode, first line
	EntryPoint *uc1;
	//Uncomplicated episode, second line
	EntryPoint *uc2;
	//Severe episode
	EntryPoint *sev;
	//Non-malaria fever
	EntryPoint *nmf;
	//List of all possible CM decisions
	Decisions *decisions;
	//minAge of humnans belonging to this CM
	double minimumAgeYears;
	//maxAge of humnans belonging to this CM
	double maximumAgeYears;
	int x;

public:
	CaseManagement(DOMNode * node) :x(0){
		createNode(this,node);
	}
	
	~CaseManagement(void) {
		delete uc1;
		delete uc2;
		delete sev;
		delete nmf;
		delete decisions;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){

		minimumAgeYears = 0.0;
		maximumAgeYears = 100.0;
		if (Converter::contains(s_MAX_AGE_YRS,map)){
			maximumAgeYears = Converter::parseDouble(s_MAX_AGE_YRS,map);
		}
		if (Converter::contains(s_MIN_AGE_YRS,map)){
			minimumAgeYears = Converter::parseDouble(s_MIN_AGE_YRS,map);
		}
	}

	void addChild(DOMNode * child){
		
		if (Converter::equals(child,s_UC1)){
			uc1 =  new EntryPoint(child);			
		}	
		if (Converter::equals(child,s_UC2)){
			uc2 =  new EntryPoint(child);			
		}
		if (Converter::equals(child,s_SEV)){
			sev =  new EntryPoint(child);			
		}
		if (Converter::equals(child,s_NMF)){
			nmf =  new EntryPoint(child);			
		}
		if (Converter::equals(child,s_DECISIONS)){
			decisions =  new Decisions(child);			
		}
	}

	EntryPoint* getUncomplicatedFirst(){
		return uc1;
	}
	EntryPoint* getUncomplicatedSecond(){
		return uc2;
	}
	EntryPoint* getSevere(){
		return sev;
	}
	EntryPoint* getNMF(){
		return nmf;
	}
	Decisions* getDecisions(){
		return decisions;
	}

	double getMinAgeYrs(){
		return minimumAgeYears;
	}

	double getMaxAgeYrs(){
		return maximumAgeYears;
	}

#ifdef _LOG
	void debug(){
		cerr << "<Case Management " <<
			"\t>\n";
	}
#endif

};

#endif
