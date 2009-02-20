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



#ifndef GENO_TYPE_H
#define GENO_TYPE_H

#include "Names.h"
#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class GenoType : public MalariaNode
{
private:		
	string name;
	double freq;
	double acr;
	int proph;
	int tolperiod;
	double atten;

public:
	GenoType(DOMNode * node){
		createNode(this,node);
	}
	
	~GenoType(void) {}
	

	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		name = Converter::getValue(s_NAME,map);
		freq = Converter::parseDouble(s_FREQ,map);
		acr = Converter::parseDouble(s_ACR,map);
		proph = Converter::parseInt(s_PROPH,map);
		tolperiod = Converter::parseInt(s_TOLPERIOD,map);
		atten = Converter::parseDouble(s_ATTEN,map);


	}

	void addChild(DOMNode * child){	
	}

	double getFreq(){
		return freq;
	}

	double getACR(){
		return acr;
	}

	int getProph(){
		return proph;
	}
	int getTolperiod(){
		return tolperiod;
	}

	double getAtten(){
		return atten;
	}


#ifdef _LOG
	void debug(){
		cerr << "<GenoType " <<
			"\t name \t" << name <<
			"\t freq \t" << freq <<
			"\t acr \t" << acr <<
			"\t proph \t" << proph <<
			"\t tolperiod \t" << tolperiod <<
			"\t atten \t" << atten <<
			"\t>\n";
	}
#endif

};

#endif

