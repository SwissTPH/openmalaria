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



#ifndef STRING_NODE_CPP
#define STRING_NODE_CPP

#include "Names.h"
#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

//This class just reads, the text within the element as a double value. 
//This class does not delete the char, you will have to delete it
class StringNode : public MalariaNode
{
private:
	char * value;
public:
	StringNode(DOMNode * node){
		createNode(this,node);
	}
	
	
	~StringNode(void) {	
		//Do not delete the name, because caller will do it
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		for (int i = 0, length = list->getLength(); i < length;  i++){
			DOMNode * child = list->item(i);
			if (child->getNodeType() == DOMNode::TEXT_NODE)
				value = XMLString::transcode(child->getNodeValue());
		}
	}

	void addChild(DOMNode * child){	
		
	}

	char * getValue(){
		return value;
	}

#ifdef _LOG
	void debug(){}
#endif

};

#endif
