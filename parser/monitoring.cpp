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


#ifndef MONITORING_CPP
#define MONITORING_CPP

#include "Names.h"

#include "converter.cpp"
#include "surveys.cpp"
#include "ageGroup.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Monitoring : public MalariaNode{
private:
	Surveys* surveys;
	AgeGroup* ageGroup;
public:
	Monitoring(DOMNode * node){
		createNode(this,node);
	};

	~Monitoring(void) {
		delete surveys;
		delete ageGroup;
	};
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
	};

	void addChild(DOMNode * child){
		
		if (Converter::equals(child, s_SURVEYS)){
			surveys =  new Surveys(child);			
		}		
		
		if (Converter::equals(child, s_AGE_GROUP)){
			ageGroup =  new AgeGroup(child);			
		}		
	}

	bool isSurvey(int time){
		return surveys->isSurvey(time);
	}

	Surveys * getSurveys(){
		return surveys;
	}

	AgeGroup * getAgeGroup(){
		return ageGroup;
	}

#ifdef _LOG
	void debug(){
		cerr << "<monitoring " <<
			"\t>\n";
	}
#endif

};


#endif
