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

#ifndef IPT_DESCRIPTION_CPP
#define IPT_DESCRIPTION_CPP

#include "genoType.cpp"

#include "Names.h"
#include "value.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class IPTDescription : public MalariaNode
{
private:	
	Value * halfLifeYrs;
	double iptiEffect;

	int numGenoTypes;
	GenoType **genotypes;
	int x;

public:
	IPTDescription(DOMNode * node) :x(0){
		createNode(this,node);
	}
	
	~IPTDescription(void) {
		delete halfLifeYrs;

		for (int i =0; i < numGenoTypes; i++)
			delete genotypes[i];
		delete [] genotypes;

	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		//Get the iptiEffect: the value is a double
		iptiEffect =  Converter::parseDouble(s_IPTIEFFECT,map);				
		//Check the number of Genotypes

		numGenoTypes = 0;
		for (int i = 0,length=list->getLength() ; i < length ; i++){
			DOMNode * node = list->item(i);
			if (node->getNodeType() == DOMNode::ELEMENT_NODE){
				if (Converter::equals(node,s_INFGENOTYPE))				
					numGenoTypes ++;
			}
		}
		//We can create the array
		genotypes = new GenoType*[numGenoTypes];

	}

	void addChild(DOMNode * child){	
		if (Converter::equals(child,s_HALF_LIFE_YRS)){			
			halfLifeYrs = new Value(s_HALF_LIFE_YRS,child);
		}

		if (Converter::equals(child,s_INFGENOTYPE)){
			genotypes[x++] = new GenoType(child);
		}		
	}

	Value * getHalfLifeYrs(){
		return halfLifeYrs;
	}

	int getnumGenoTypes(){
		return numGenoTypes;
	}

	double getIptiEffect(){
		return iptiEffect;
	}

	GenoType * getGenoType(int index){
		return genotypes[index];
	}


#ifdef _LOG
	void debug(){		
		cerr << "<IPTIDescription " <<
			"\tITPIEffect " << iptiEffect <<
			"\t>\n";
}
#endif

};

#endif
