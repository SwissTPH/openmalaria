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



#ifndef MEDICATE_H
#define MEDICATE_H

#include "Names.h"
#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class Medicate : public MalariaNode{
private:	
	double qty;
	int time;
	char * name;
	
public:
	Medicate(DOMNode * node){
		createNode(this,node);
	}
	
	~Medicate(void) {}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		qty = Converter::parseDouble(s_QTY,map);
		time = Converter::parseInt(s_TIME,map);
		if (Converter::contains(s_NAME,map)){
			name = Converter::getValue(s_NAME,map);
	}
	}

	void addChild(DOMNode * child){	
	}

	double getQTY(){
		return qty;
	}

	int getTime(){
		return time;
	}
	
	char* getName() {
        return name;
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

