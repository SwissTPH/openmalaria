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

#ifndef ITEM_CPP
#define ITEM_CPP

#include "Names.h"
#include "converter.cpp"

#include "value.cpp"
#include <iostream>
using namespace std;

class Item : public MalariaNode{
private:
  double maxAgeYrs;
  double value;

public:
  Item(DOMNode * node){      
    createNode(this,node);
  }
  
  ~Item(void){
  }
  
  void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
    if (Converter::contains(s_VALUE,map)){
       value = Converter::parseDouble(s_VALUE,map);
    }else{
      value = Converter::parseDouble(s_BEST,map);
    }
    maxAgeYrs = Converter::parseDouble(s_MAX_AGE_YRS,map);
  }

  void addChild(DOMNode * child){  
  }

  double getValue(){
    return value;
  }

  double getMaxAgeYrs(){
    return maxAgeYrs;
  }
  
#ifdef _LOG
  void debug(){
    cerr << "<item " <<
      "\tMaxAgeYrs " << maxAgeYrs <<      		
      "\t>\n";    
  }
#endif

};

#endif
