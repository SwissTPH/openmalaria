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



#ifndef ENTODATA_CPP
#define ENTODATA_CPP

#include "converter.cpp"
#include "EIRDaily.cpp"
#include "anopheles.cpp"

#include "MalariaNode.h"
#include "Constant.h"
#include <iostream>
using namespace std;

#include "Names.h"

enum InputType {EIR,FOURIER};

class EntoData : public MalariaNode{
private:	
  InputType type;
  char * name;
  int numEir;
  EIRDaily ** dailies;
  int x;
  Anopheles * anopheles;
    
public:
  EntoData(DOMNode * node) : x(0){
    createNode(this,node);
  };
  ~EntoData(void){
    delete [] name;
    for (int i =0; i < numEir; i++){
      delete dailies[i];
    }
      delete [] dailies;
  };
	
  void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
    name = Converter::getValue(s_NAME, map);
    //Compare the type: EIR or VC into the XMLTranscode: it is faster
    XMLCh * xml_eir = XMLString::transcode(s_EIR); 
    XMLCh * xml_inputtype = XMLString::transcode(s_INPUTTYPE);
    type = XMLString::equals(
             map->getNamedItem(xml_inputtype)->getNodeValue(),xml_eir) ? EIR : FOURIER;
    XMLString::release(&xml_eir);
    XMLString::release(&xml_inputtype);
    numEir = 0;
    if (type==EIR){
      for (int i = 0, length = list->getLength(); i < length ; i++){
        if (list->item(i)->getNodeType() == DOMNode::ELEMENT_NODE){
          numEir++;
        }
      }
      dailies = new EIRDaily* [numEir];		
    }
  }

  void addChild(DOMNode * child){
    if (type==EIR){	
      dailies[x++] = new EIRDaily(child);
    }
    else{
      anopheles= new Anopheles(child);
    }
  }

  double getEIRDaily(int index){
    return (index >= 0 && index < numEir)? dailies[index]->getValue() : MISSING_VALUE;
  }

  Anopheles* getAnopheles(string name){
  //TODO: return Anopheles by name
  return anopheles;
  }

#ifdef _LOG
	void debug(){
		cerr << "<entoData" ;
		if (name != NULL)	
			cerr << "\tname " << name;
		cerr << "\tinputType " << type <<
			"\t>\n";
	}
#endif

};

#endif
