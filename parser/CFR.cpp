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

#ifndef CFR_CPP
#define CFR_CPP

#include "CFRGroup.cpp"

#include "Names.h"
#include "converter.cpp"


#include "MalariaNode.h"
#include <iostream>

using namespace std;

class CFR : public MalariaNode{
private:	
	int numGroups;
	CFRGroup **groups;
	int x;
public:
	CFR(DOMNode * node) :x(0){
		createNode(this,node);
	}
	
	~CFR(void) {
		for (int i =0; i < numGroups; i++)
			delete groups[i];
		delete [] groups;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		numGroups = 0;
		for (int i = 0,length=list->getLength() ; i < length ; i++){
			DOMNode * node = list->item(i);
			if (node->getNodeType() == DOMNode::ELEMENT_NODE){
				if (Converter::equals(node,s_GROUP))				
					numGroups ++;
			}
		}
		groups = new CFRGroup*[numGroups];
	}

	void addChild(DOMNode * child){	
	  if (Converter::equals(child,s_GROUP))
	    groups[x++] = new CFRGroup(child);
	  }

	//return the number of agegroups for which we have CFR
	int getNumGroups(){
		return numGroups;
	}

	CFRGroup * getGroup(int index){
		return groups[index];
	}

#ifdef _LOG
	void debug(){}
#endif

};


#endif

