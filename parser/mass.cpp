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



#ifndef MASS_CPP
#define MASS_CPP

#include "Names.h"
#include "Constant.h"

#include "converter.cpp"
#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Mass : public MalariaNode{
private:
  const char * name;
  double maxAge;
  double minAge;
  double coverage;
  bool b_maxAge;
  bool b_minAge;
  bool b_coverage;
public:
	Mass(DOMNode * node, const char *objectName) : name(objectName){
                XMLString::transcode(s_MAXIMUM_AGE_YEARS);
		createNode(this,node);		
	};

	~Mass(void) {		
	};
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
	  //Sets the default values;
	  minAge = MIN_AGE;
	  maxAge = MAX_AGE;
	  coverage = COMPLIANCE;		
	  //TODO: This is not very efficient.
          XMLCh * xml_min_age=XMLString::transcode(s_MIN_AGE);
	  XMLCh * xml_max_age=XMLString::transcode(s_MAX_AGE);
	  XMLCh * xml_coverage=XMLString::transcode(s_COVERAGE);
	  if ((b_minAge = map->getNamedItem(xml_min_age)) != false){
	    minAge = Converter::parseDouble(s_MIN_AGE,map);
	  }
		
	  if ((b_maxAge = map->getNamedItem(xml_max_age)) != false){
	    maxAge = Converter::parseDouble(s_MAX_AGE,map);
	  }		
	
	  if ((b_coverage = map->getNamedItem(xml_coverage)) != false){
	    coverage = Converter::parseDouble(s_COVERAGE,map);
	  }				
          XMLString::release(&xml_min_age);
	  XMLString::release(&xml_max_age);
	  XMLString::release(&xml_coverage);
	}

	void addChild(DOMNode * child){	
		
	}

	bool isMaxAge(){
		return b_maxAge;
	}

	bool isMinAge(){
		return b_minAge;
	}

	bool isCoverage(){
		return b_coverage;
	}

	double getMinAge(){
		return minAge;
	}

	double getMaxAge(){
		return maxAge;
	}

	double getCoverage(){
		return coverage;
	}

#ifdef _LOG
	void debug(){
		cerr << "<" << name ;
		if (b_minAge)
			cerr << "\tminAge " << minAge ;
		if (b_maxAge)
			cerr << "\tmaxAge " << maxAge ;
		if (b_coverage)
			cerr << "\tcoverage " << coverage;
		cerr << "\t>\n";
	}
#endif

};

#endif
