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

#ifndef INTERVENTION_CPP
#define INTERVENTION_CPP

#include "Names.h"

#include "converter.cpp"
#include "entoData.cpp"
#include "healthSystem.cpp"
#include "mass.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;



class Intervention : public MalariaNode
{
private:	
	int time;
	Mass * vaccinate;
	Mass * MDA ;
	EntoData * changeEIR;
	HealthSystem * changeHS;
	//TODO Not used
	Mass * ITN ;
	Mass * IPTi ;
	bool b_vaccinate;
	bool b_MDA;
	bool b_changeEIR;
	bool b_changeHS;
	bool b_IRS;
	bool b_ITN;
	bool b_IPTi;

public:
	Intervention(DOMNode * node)
	{
		createNode(this,node);
	}

	~Intervention(void) {
		if (b_vaccinate) delete vaccinate;
		if (b_MDA) delete MDA;		
		if (b_ITN) delete ITN;
		if (b_IPTi) delete IPTi;

		//Don't need to delete changeEIR or changeHS. This happens in InputDataBoinc.
		//TODO: check this
		//if (b_changeEIR) delete changeEIR;
		//if (b_changeHS) delete changeHS;
	};
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){			
		b_vaccinate = false;
		b_MDA = false;
		b_IRS = false;
		b_changeEIR = false;
		b_changeHS = false;
		b_ITN = false;
		b_IPTi = false;
		time = Converter::parseInt(s_TIME,map);
	}

	void addChild(DOMNode * child){			
		
		if (Converter::equals(child,s_VACCINATE)){
			vaccinate = new Mass(child,s_VACCINATE);
			b_vaccinate = true;
		}		
		if (Converter::equals(child,s_MDA)){
			MDA = new Mass(child,s_MDA);
			b_MDA = true;
		}
		if (Converter::equals(child,s_ITN)){
			ITN = new Mass(child,s_ITN);
			b_ITN = true;
		}
		if (Converter::equals(child,s_CHANGE_EIR)){
			changeEIR = new EntoData(child);
			b_changeEIR = true;
		}
		if (Converter::equals(child,s_CHANGE_HS)){
			changeHS = new HealthSystem(child);
			b_changeHS = true;
		}
		if (Converter::equals(child,s_IRS)){			
			b_IRS = true;
		}
		if (Converter::equals(child,s_IPTI)){		
			IPTi = new Mass(child,s_IPTI);
			b_IPTi = true;
		}

	}

	int getTime(){
		return time;
	}

	bool isVaccine(){
		return b_vaccinate;
	}

	bool isIRS(){
		return b_IRS;
	}

	bool isMDA(){
		return b_MDA;
	}

	bool isITN(){
		return b_ITN;
	}

	bool isChangeEIR(){
		return b_changeEIR;
	}

	bool isChangeHS(){
		return b_changeHS;
	}

	bool isIPTi(){
		return b_IPTi;
	}

	Mass * getMDA(){
		return MDA;
	}

	Mass * getITN(){
		return ITN;
	}

	Mass * getVaccine(){
		return vaccinate;
	}

	EntoData * getChangeEIR(){
		return changeEIR;
	}

	HealthSystem * getChangeHS(){
		return changeHS;
	}

	Mass * getIPTi(){
		return IPTi;
	}

#ifdef _LOG
	void debug(){
		cerr << "<Intervention " << 			
			"\ttime " << time <<
			"\t>\n";
	}
#endif

};

#endif
