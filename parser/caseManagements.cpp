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



#ifndef CASEMANAGEMENTS_H
#define CASEMANAGEMENTS_H

#include "Names.h"
#include "caseManagement.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class CaseManagements : public MalariaNode{
private:	
	CaseManagement **caseManagements;
	int numCaseManagements;
	//CaseManagement* cm;
	int x;

public:
	CaseManagements(DOMNode * node):x(0){
		createNode(this,node);
	}
	
	~CaseManagements(void){
		for (int i =0; i < numCaseManagements; i++)
			delete caseManagements[i];
		delete [] caseManagements;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		numCaseManagements = 0;
		for (int i = 0,length=list->getLength() ; i < length ; i++){
			DOMNode * node = list->item(i);
			if (node->getNodeType() == DOMNode::ELEMENT_NODE){
				if (Converter::equals(node,s_CASE_MANAGEMENT))				
					numCaseManagements ++;
			}
		}
		caseManagements = new CaseManagement*[numCaseManagements];
	}

	void addChild(DOMNode * child){	
		if (Converter::equals(child,s_CASE_MANAGEMENT))
			caseManagements[x++] = new CaseManagement(child);
	}

	//Note the difference in the index between here and addChild() above. All 1-based to 0-based
	//index translation happens in inputData.cpp
	CaseManagement* getCaseManagementbyID(int id){
		return caseManagements[id];
	}

	CaseManagement* getCaseManagementbyAge(double Age){
		for (int i =0; i < numCaseManagements; i++){
			if (caseManagements[i]->getMinAgeYrs() <= Age && Age <= caseManagements[i]->getMaxAgeYrs())
				return caseManagements[i];
		}
		//default TODO: check if this is a good thing to return
		cerr << "No Case Management for Age Class found";
		return caseManagements[0];
	}



#ifdef _LOG
	void debug(){
	}
#endif

};

#endif

