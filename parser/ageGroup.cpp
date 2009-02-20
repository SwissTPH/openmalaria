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


#ifndef AGE_GROUP_CPP
#define AGE_GROUP_CPP
/*test*/
#include "Names.h"

#include "group.cpp"

#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class AgeGroup : public MalariaNode{
private:

	//The lowerbound of the age
	double lowerBound;	
	//The different groups 
	Group ** groups;
	//The number of groups
	int numGroups;
	int x;

public:

	AgeGroup(DOMNode * node):x(0){
		createNode(this,node);
	}

	~AgeGroup(void) {
		//We have to delete each group
		for (int i = 0; i < numGroups; i++)
			delete groups[i];
		//We can now delete the array
		delete [] groups;
	}

	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		//Get the lower bound: the value is a double
		lowerBound =  Converter::parseDouble(s_LOWER_BOUND,map);				
		//Check the number of groups
		numGroups = 0;
		for (int i =0, length = list->getLength(); i < length ; i++){
			DOMNode * child = list->item(i);
			//We have to check if it is really an element
			if (child->getNodeType() == DOMNode::ELEMENT_NODE)
				numGroups++;
		}
		//We can create the array
		groups = new Group*[numGroups];
	};

	//Adds a new group: the method is called by its superclass
	//There can be only one element group.
	void addChild(DOMNode * child){				
		groups[x++] = new Group(child);
	}

	//Return the lowerbound of the ageGroup
	double getLowerBound(){
		return lowerBound;
	}

	//return the number of groups
	int getNumGroups(){
		return numGroups;
	}

	//return the (index)th group
	Group * getGroup(int index){
		return groups[index];
	}

#ifdef _LOG
	void debug(){
		cerr << "<ageGroup " <<  
			"\tlowerbound " << lowerBound << 
			"\tnumGroups " << numGroups <<  
			"\t>\n";
	}
#endif
};

#endif
