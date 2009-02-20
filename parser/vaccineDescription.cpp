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



#ifndef VACCINE_DESCRIPTION_CPP
#define VACCINE_DESCRIPTION_CPP

#include "converter.cpp"
#include "value.cpp"
#include "Names.h"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class VaccineDescription : public MalariaNode
{
private:	
	int vaccineType;
	Value * halfLifeYears;
	Value * efficacyB;
	int numInitialEfficacy;
	Value ** initialEfficacies;
	int current;
public:
	VaccineDescription(DOMNode * node): current(0){
		createNode(this,node);		
	}
	
	~VaccineDescription(void) {
		delete halfLifeYears;
		delete efficacyB;
		for (int i = 0; i < numInitialEfficacy; i++)
			delete initialEfficacies[i];
		delete [] initialEfficacies;				
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		vaccineType = Converter::parseInt(s_VACCINE_TYPE,map);
		//Count the number of initialEfficacy
		numInitialEfficacy = 0;
		for (int i =0, length = list->getLength(); i < length ;i++){
			DOMNode * node = list->item(i);
			if (Converter::equals(node,s_INITIAL_EFFICACY)){
					numInitialEfficacy++;
				}
			initialEfficacies = new Value*[numInitialEfficacy];
		}
	}

	void addChild(DOMNode * child){	
		if (Converter::equals(child,s_HALF_LIFE_YRS)){
			halfLifeYears = new Value(s_HALF_LIFE_YRS,child);
		}
                if (Converter::equals(child,s_EFFICACY_B)){
			efficacyB = new Value(s_EFFICACY_B,child);
		}
                if (Converter::equals(child,s_INITIAL_EFFICACY)){
			initialEfficacies[current++] = new Value(s_INITIAL_EFFICACY,child);
		}
	}

	int getVaccineType(){
		return vaccineType;
	}

	Value * getInitialEfficacy(int index){
		return initialEfficacies[min(index,numInitialEfficacy-1)];
	}

	int getNumInitialEfficacy(){
		return numInitialEfficacy;
	}

	Value * getHalfLifeYears(){
		return halfLifeYears;
	}

	Value * getEfficacyB(){
		return efficacyB;
	}
	
#ifdef _LOG
	void debug(){
		cerr << "<vaccineDescription " <<
			"\tvaccineType " << vaccineType <<
			"\t>\n";
	}
#endif

};

#endif

