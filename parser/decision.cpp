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



#ifndef DECISION_H
#define DECISION_H

#include "Names.h"
#include "converter.cpp"
#include "medicate.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class Decision : public MalariaNode{
private:	
	int id;
	Medicate **medicates;
	int numMedicates;
	int x;
	
public:
	Decision(DOMNode * node):x(0){
		createNode(this,node);
	}
	
	~Decision(void){
			for (int i =0; i < numMedicates; i++)
			delete medicates[i];
		delete [] medicates;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		id = Converter::parseInt(s_ID,map);
		numMedicates = 0;
		for (int i = 0,length=list->getLength() ; i < length ; i++){
			DOMNode * node = list->item(i);
			if (node->getNodeType() == DOMNode::ELEMENT_NODE){
				if (Converter::equals(node,s_MEDICATE))				
					numMedicates ++;
			}
		}
		medicates = new Medicate*[numMedicates];
	}

	void addChild(DOMNode * child){	
		if (Converter::equals(child,s_MEDICATE))
			medicates[x++] = new Medicate(child);
	}

	int getID(){
		return id;
	}

	int getNumMedicates(){
		return numMedicates;
	}

	Medicate* getMedicate(int index){
		return medicates[index];
	}

#ifdef _LOG
	void debug(){
	}
#endif

};

#endif

