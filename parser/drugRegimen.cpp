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



#ifndef DRUG_REGIMEN_CPP
#define DRUG_REGIMEN_CPP

#include "Names.h"
#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class DrugRegimen : public MalariaNode
{
private:	
	char * firstLine;
	char * secondLine;
	char * inpatient;

public:
	DrugRegimen(DOMNode * node){
		createNode(this,node);
	}

	~DrugRegimen(void) {
		delete firstLine;
		delete secondLine;
		delete inpatient;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		firstLine = Converter::getValue(s_FIRST_LINE,map);
		secondLine = Converter::getValue(s_SECOND_LINE,map);
		inpatient = Converter::getValue(s_INPATIENT,map);
	}

	void addChild(DOMNode * child){	
	}

	char * getFirstLine(){
		return firstLine;
	}

	char * getSecondLine(){
		return secondLine;
	}

	char * getInpatient(){
		return inpatient;
	}


#ifdef _LOG
	void debug(){
		cerr << "<DrugRegimen " <<
			"\tfirstLine " << firstLine <<
			"\tsecondLine " << secondLine <<
			"\tinpatient " << inpatient<<
			"\t>\n";
	}
#endif

};

#endif
