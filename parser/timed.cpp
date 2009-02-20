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



#ifndef TIMED_CPP
#define TIMED_CPP

#include "intervention.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;



class Timed : public MalariaNode
{
private:	
	Intervention ** interventions;
	int numGroups;
	bool sort;
	int x;
public:
	Timed(DOMNode * node):x(0){
		createNode(this,node);
	}

	~Timed(void) {	
		for (int i = 0; i < numGroups; i++)
			delete interventions[i];
		delete [] interventions;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		//Check the number of nodes
		sort =false;
		numGroups = 0;
		for (int i =0, length = list->getLength(); i < length ; i++){
			DOMNode * child = list->item(i);
			if (child->getNodeType() == DOMNode::ELEMENT_NODE)
				numGroups++;
		}
		interventions = new Intervention*[numGroups];
	}

	void addChild(DOMNode * child){					
		interventions[x++] = new Intervention(child);
		if (x > 1 && 
			interventions[x-2]->getTime() >= interventions[x-1]->getTime()) 
			sort = false;
	}

	Intervention * getIntervention(int index){
		//It is possible that we need to check the size...
		return interventions[index];
	}

	Intervention * getInterventionByTime(int time){
		if (sort){
			//Dichotomic search
			int minIndex = 0;
			int maxIndex = numGroups-1;
			if (interventions[minIndex]->getTime() == time) 
				return interventions[minIndex];
			if (interventions[maxIndex]->getTime() == time) 
				return interventions[maxIndex];
			//We can decrement maxIndex and increment minIndex
			//Because we know that there are not at the border
			while (minIndex+1 < maxIndex){
				int middle = (maxIndex+minIndex)/2;
				int value = interventions[middle]->getTime();
				if (value == time) return interventions[middle];
				if (value < time)
					maxIndex = middle;
				else 
					minIndex = middle;
			}
			return NULL;
		} else {
			//Brute force search
			for (int i = 0; i < numGroups; i++){
				if (time == interventions[i]->getTime())
					return interventions[i];
			}
			return NULL;
		}
	}

	int getNumInterventions(){
		return numGroups;
	}

#ifdef _LOG
	void debug(){
		cerr << "<timed " << "\t" <<
			"\tsort " << sort <<
			"\tnumInterventions " << numGroups << "\t" <<
			">\n";
	}
#endif

};

#endif
