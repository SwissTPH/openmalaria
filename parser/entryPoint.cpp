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



#ifndef ENTRY_POINT_H
#define ENTRY_POINT_H

#include "Names.h"
#include "converter.cpp"
#include "endPoint.cpp"

#include "GSLWrapper.h"
#include "MalariaNode.h"
#include <iostream>
using namespace std;

class EntryPoint : public MalariaNode{
private:	
	int numEndPoint;
	EndPoint **endPoints;

public:
	EntryPoint(DOMNode * node){
		createNode(this,node);
	}
	
	~EntryPoint(void) {
		for (int i =0; i < numEndPoint; i++)
			delete endPoints[i];
		delete [] endPoints;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		numEndPoint = 0;
		for (int i = 0,length=list->getLength() ; i < length ; i++){
			DOMNode * node = list->item(i);
			if (node->getNodeType() == DOMNode::ELEMENT_NODE){
				if (Converter::equals(node,s_ENDPOINT))				
					numEndPoint++;
			}
		}
		endPoints = new EndPoint*[numEndPoint];
	}

	void addChild(DOMNode * child){	
		if (Converter::equals(child,s_ENDPOINT)){
			EndPoint* ep = new EndPoint(child);
			endPoints[ep->getDecision()-1] = ep;
		}
	}

	int getDecisionID(){
		double ranNum=W_UNIFORM();
		double cumulP=0;
		for (int i =0; i < numEndPoint; i++){
			cumulP+=endPoints[i]->getP();
			if (cumulP>=ranNum){
				return endPoints[i]->getDecision();
			}
		}
		cerr << "No oc id:" << cumulP << ranNum << "\n";
		return endPoints[numEndPoint-1]->getDecision();
	}

#ifdef _LOG
	void debug(){
		cerr << "<group " <<
			"\tcfr " << cfr <<
			"\tlowerbound " << lowerbound <<
			"\t>\n";
	}
#endif

};

#endif

