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


#ifndef PARAMETERS_CPP
#define PARAMETERS_CPP

#include "parameter.cpp"
#include "converter.cpp"

#include "Names.h"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class Parameters : public MalariaNode
{
private:	
	double latentp;
	double delta;
	double nspore;
	int interval;
	int iSeed;

	int numParameters;
	int x;

	Parameter ** parameters;
public:
	Parameters(DOMNode * node) : x(0){
		createNode(this,node);
	}
	
	~Parameters(void) {
		for (int i = 0; i < numParameters; i++)
			delete parameters[i];
		delete [] parameters;
	}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		latentp = Converter::parseDouble(s_LATENT_P,map);
		delta = Converter::parseDouble(s_DELTA,map);
		nspore = Converter::parseDouble(s_NSPORE,map);
		interval = Converter::parseInt(s_INTERVAL,map);
		iSeed = Converter::parseInt(s_I_SEED,map);
		numParameters = 0;
		for (int i =0, length = list->getLength() ; i < length ; i++){
			DOMNode * node = list->item(i);
			if (node->getNodeType() == DOMNode::ELEMENT_NODE){
				if (Converter::equals(node,s_PARAMETER))
					numParameters++;
			}
		}
		parameters = new Parameter*[numParameters];
	}

	void addChild(DOMNode * child){	
		parameters[x++] = new Parameter(child);
	}

	double getLatentp(){
		return latentp;
	}

	double getDelta(){
		return delta;
	}

	double getNspore(){
		return nspore;
	}

	int getISeed(){
		return iSeed;
	}

	int getInterval(){
		return interval;
	}

	Parameter * getParameter(int index){
		return parameters[index];
	}

	int getNumParameters(){
		return numParameters;
	}

#ifdef _LOG
	void debug(){
		cerr << "<Parameters " <<
			"\numParams " << numParameters <<
			"\tlatentp " << latentp <<
			"\tdelta " << delta <<
			"\tinterval " << interval <<
			"\tiSeed " << iSeed <<
			"\t>\n";
	}
#endif

};

#endif

