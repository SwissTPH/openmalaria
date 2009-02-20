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

#ifndef ITN_DESCRIPTION_CPP
#define ITN_DESCRIPTION_CPP

#include "Names.h"
#include "value.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class ITNDescription : public MalariaNode
{
private:	
	Value * halfLifeYrs;
	Value * pu0;
	Value * pu1;
	Value * sporogonyGonotrophy;
public:
	ITNDescription(DOMNode * node){
		createNode(this,node);
	}
	
	~ITNDescription(void) {
		delete halfLifeYrs;
		delete pu0;
		delete pu1;
		delete sporogonyGonotrophy;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	

	}

	void addChild(DOMNode * child){	
	  if (Converter::equals(child,s_HALF_LIFE_YRS)){			
	    halfLifeYrs = new Value(s_HALF_LIFE_YRS,child);
	  }		
		
	  if (Converter::equals(child,s_PU_0)){			
	    pu0 = new Value(s_PU_0,child);
	  }		
	
	  if (Converter::equals(child,s_PU_1)){			
	    pu1 = new Value(s_PU_1,child);
	  }		
				
	  if (Converter::equals(child,s_SPOROGONY_GONOTROPHY)){			
	    sporogonyGonotrophy = new Value(s_SPOROGONY_GONOTROPHY,child);
	  }				
	}

	Value * getPu0(){
		return pu0;
	}

	Value * getPu1(){
		return pu1;
	}

	Value * getSporogonyGonotrophy(){
		return sporogonyGonotrophy;
	}

	Value * getHalfLifeYrs(){
		return halfLifeYrs;
	}

#ifdef _LOG
	void debug(){}
#endif

};

#endif
