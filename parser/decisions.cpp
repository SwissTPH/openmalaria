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



#ifndef DECISIONS_H
#define DECISIONS_H

#include "Names.h"
#include "decision.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class Decisions : public MalariaNode{
private:	
	Decision **decisions;
	int numDecisions;

public:
	Decisions(DOMNode * node){
		createNode(this,node);
	}
	
	~Decisions(void){
		for (int i =0; i < numDecisions; i++)
			delete decisions[i];
		delete [] decisions;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		numDecisions = 0;
		for (int i = 0,length=list->getLength() ; i < length ; i++){
			DOMNode * node = list->item(i);
			if (node->getNodeType() == DOMNode::ELEMENT_NODE){
				if (Converter::equals(node,s_DECISION))				
					numDecisions ++;
			}
		}
		decisions = new Decision*[numDecisions];
	}

	void addChild(DOMNode * child){	
		if (Converter::equals(child,s_DECISION)){
			Decision* dec = new Decision(child);
			decisions[dec->getID()-1] = dec;
		}
	}

	//Note the difference in the index between here and addChild() above. All 1-based to 0-based
	//index translation happens in inputData.cpp
	Decision* getDecision(int id){
		return decisions[id];
	}

#ifdef _LOG
	void debug(){
	}
#endif

};

#endif

