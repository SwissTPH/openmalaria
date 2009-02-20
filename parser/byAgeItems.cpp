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

#ifndef BY_AGE_ITEMS_CPP
#define BY_AGE_ITEMS_CPP

#include "item.cpp"
#include "Constant.h"

#include "Names.h"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class ByAgeItems : public MalariaNode{
private:	
	int numItems;
	Item ** items;
	int x;

public:
	ByAgeItems(DOMNode * node):x(0){		
		createNode(this,node);
	}
	
	~ByAgeItems(void) {}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		numItems = 0;
		for (int i =0,length = list->getLength(); i < length ; i++){
			DOMNode * node = list->item(i);
			if (node->getNodeType() ==  DOMNode::ELEMENT_NODE){				
				if (Converter::equals(node,s_ITEM))					
					numItems++;
			}	
		}
		items = new Item*[numItems];
	}

	void addChild(DOMNode * child){			
	if (Converter::equals(child,s_ITEM))
	  items[x++] = new Item(child);
	}

	double getByAge(double age){
	  for (int i = 0; i < numItems; i++){
	    if (items[i]->getMaxAgeYrs() > age)
	      return items[i]->getValue();
	    }
	  return MISSING_VALUE;
	}

#ifdef _LOG
	void debug(){}
#endif

};

#endif
