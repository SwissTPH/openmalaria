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


#ifndef DOCUMENT_PARAMETER_CPP
#define DOCUMENT_PARAMETER_CPP
#include "parameters.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class DocumentParameter : public MalariaNode
{
private:
	Parameters*  parameters;
public:
	DocumentParameter(DOMNode * node){
		createNode(this,node);
	}

	~DocumentParameter(void) {
		delete parameters;
	}

	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
	//Nothing to do
	};

	void addChild(DOMNode * child){
		if (child->getNodeType() == DOMNode::ELEMENT_NODE)
			parameters =  new Parameters(child);		
	}

	Parameters * getParameters() const {
		return parameters;
	}

#ifdef _LOG
	void debug(){}
#endif

};

#endif
