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



#ifndef TREATMENT_DETAILS_CPP
#define TREATMENT_DETAILS_CPP

#include "Names.h"
#include "value.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class TreatmentDetails : public MalariaNode
{
private:
	const char * name;
	bool b_CQ;
	bool b_SP;
	bool b_AQ;
	bool b_SPAQ;
	bool b_ACT;
	bool b_QN;	
	Value * CQ;
	Value * SP;
	Value * AQ;
	Value * SPAQ;
	Value * ACT;
	Value * QN;
	Value * selfTreatment;
public:
	TreatmentDetails(const char * currentName, DOMNode * node) : name(currentName){					
		createNode(this,node);
	};

	~TreatmentDetails(void) {
		if (b_CQ) delete CQ;
		if (b_SP) delete SP;
		if (b_AQ) delete AQ;
		if (b_SPAQ) delete SPAQ;
		if (b_ACT) delete ACT;
		if (b_QN) delete QN;
		delete selfTreatment;		
	};
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		b_CQ = false;
		b_SP = false;
		b_AQ = false;
		b_SPAQ = false;
		b_ACT = false;
		b_QN = false;
	}

	void addChild(DOMNode * child){	
		if (Converter::equals(child,s_CQ)){
			CQ = new Value(s_CQ,child);
			b_CQ = true;
		}

		if (Converter::equals(child,s_SP)){
			SP = new Value(s_SP,child);
			b_SP = true;
		}

		if (Converter::equals(child,s_AQ)){
			AQ = new Value(s_AQ,child);
			b_AQ = true;
		}

		if (Converter::equals(child,s_SPAQ)){
			SPAQ = new Value(s_SPAQ,child);
			b_SPAQ = true;
		}

		if (Converter::equals(child,s_ACT)){
			ACT = new Value(s_ACT,child);
			b_ACT = true;
		}

		if (Converter::equals(child,s_QN)){
			QN = new Value(s_QN,child);
			b_QN = true;
		}

		if (Converter::equals(child,s_SELF_TREATMENT)){
			selfTreatment = new Value(s_SELF_TREATMENT,child);
		}
	}

	Value * getCQ(){
		return CQ;
	}

	Value * getSP(){
		return SP;
	}

	Value * getAQ(){
		return AQ;
	}

	Value * getSPAQ(){
		return SPAQ;
	}

	Value * getACT(){
		return ACT;
	}

	Value * getQN(){
		return QN;
	}

	Value * getSelfTreatment(){
		return selfTreatment;
	}

	double getACRByName(char * name){
		if (strcmp(name,s_CQ) == 0)
			return CQ->getValue();
		
		if (strcmp(name,s_SP) == 0)
			return SP->getValue();		

		if (strcmp(name,s_AQ) == 0)
			return AQ->getValue();

		if (strcmp(name,s_ACT) == 0)
			return ACT->getValue();

		if (strcmp(name,s_QN) == 0)
			return QN->getValue();

		if (strcmp(name,s_SELF_TREATMENT) == 0){
			return selfTreatment->getValue();
		}
		return -1;
		//TODO problem here
	}

	
#ifdef _LOG
	void debug(){
		cerr << "<" << name << "\t>\n";			
	}
#endif

};

#endif
