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



#ifndef GROUP_CPP
#include "Names.h"

#include "converter.cpp"
#include "Constant.h"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Group : public MalariaNode {
private:
	double upperBound;	
	double popPercent;	
	bool hasFieldUpperBound;
public:
	Group(DOMNode * node){
		createNode(this,node);
	};
	
	~Group(void) {};
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		upperBound = MISSING_VALUE;
		popPercent = MISSING_VALUE;
		hasFieldUpperBound = Converter::contains(s_UPPER_BOUND,map);
		if (hasFieldUpperBound)
			upperBound = Converter::parseDouble(s_UPPER_BOUND,map);		
		bool hasFieldPopPercent = Converter::contains(s_POP_PERCENT,map);
		if (hasFieldPopPercent)
			popPercent = Converter::parseDouble(s_POP_PERCENT,map);		
	};

	void addChild(DOMNode * child){
		//Nothing to do
	}

	double getUpperBound(){		
		return upperBound;
	}

	double getPopPercent(){
		return popPercent;
	}

#ifdef _LOG
	void debug(){
		cerr << "<group " ;
		if (hasFieldUpperBound)
			cerr << "\tupperbound " << upperBound;
		cerr << "\tpoppercent " << popPercent <<
			"\t>\n";
	}
#endif

};

#endif
