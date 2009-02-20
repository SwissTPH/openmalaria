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



#ifndef CONTINUOUS_CPP
#define CONTINUOUS_CPP

#include "Names.h"
#include "ageSpecific.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class Continuous : public MalariaNode
{
private:	
	int numVaccineDoses;
	AgeSpecific ** vaccine;

	int numItn;
	AgeSpecific ** itn;

	int numIpti;
	AgeSpecific ** ipti;

	int i ; //itn
	int v ; //vaccine
	int p ; //iptn

public:
	Continuous(DOMNode * node) :
	i(0), v(0), p(0) {
		createNode(this,node);
	}
	
	~Continuous(void) {}
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		numItn = 0;
		numVaccineDoses = 0;
		numIpti = 0;
		for (int i = 0, length = list->getLength(); i < length; i++){
			DOMNode * child = list->item(i);
			if (child->getNodeType() == DOMNode::ELEMENT_NODE){
				if (Converter::equals(child,s_VACCINE))
					numVaccineDoses++;

				if (Converter::equals(child,s_ITN))
					numItn++;

				if (Converter::equals(child,s_IPTI))
					numIpti++;

			}
		}
		vaccine = new AgeSpecific*[numVaccineDoses];
		itn = new AgeSpecific*[numItn];
		ipti = new 	AgeSpecific*[numIpti];	
	}

	void addChild(DOMNode * child){	
	  const XMLCh* nodeName = child->getNodeName();
	  XMLCh * xml_itn=XMLString::transcode(s_ITN);
          XMLCh * xml_vaccine=XMLString::transcode(s_VACCINE);
          XMLCh * xml_ipti=XMLString::transcode(s_IPTI);
		if (XMLString::equals(xml_itn,nodeName)){	
			itn[i++] = new AgeSpecific(s_ITN,child);
		}

		if (XMLString::equals(xml_vaccine,nodeName)){	
			vaccine[v++] = new AgeSpecific(s_VACCINE,child);
		}

		if (XMLString::equals(xml_ipti,nodeName)){	
			ipti[p++] = new AgeSpecific(s_IPTI,child);
		}
	  XMLString::release(&xml_itn);
	  XMLString::release(&xml_vaccine);
	  XMLString::release(&xml_ipti);
	}

	int getNumITN(){
		return numItn;
	}

	int getNumVaccineDoses(){
		return numVaccineDoses;
	}

	int getNumIPTi(){
		return numIpti;
	}


	AgeSpecific * getVaccine(int index){
		return vaccine[index];
	}

	AgeSpecific * getITN(int index){
		return itn[index];
	}

	AgeSpecific * getIPTi(int index){
		return ipti[index];
	}

#ifdef _LOG
	void debug(){}
#endif

};

#endif
