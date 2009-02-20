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


#ifndef CONVERTER_CPP
#define CONVERTER_CPP

#include <iostream>
using namespace std;
#ifdef new
#undef new
#define _REDEF_NEW
#endif
#include <xercesc/dom/DOM.hpp>
#undef new
#ifdef _REDEF_NEW
#define new DEBUG_NEW
#undef _REDEF_NEW
#endif
XERCES_CPP_NAMESPACE_USE

//In this class, each function is static inline: the goal is to accelerate the code
class Converter{
private:

public :
  //Return the atribute node that has the name "name"
  inline static DOMNode * getNode(const char * name, DOMNamedNodeMap * map){
    XMLCh * translation = XMLString::transcode(name);
    DOMNode * node = map->getNamedItem(translation);
    XMLString::release(&translation);
    return node;
  }

  //Parse the attribute name and return its integral value
  inline static int parseInt(const char * name, DOMNamedNodeMap * map){		
    return XMLString::parseInt(getNode(name,map)->getNodeValue());		
  }

  //return the value of attribute name as a string format
  inline static char * getValue(const char * name, DOMNamedNodeMap * map){
    DOMNode * node = getNode(name,map);
    if (node == NULL) return NULL;
    const XMLCh * value = node->getNodeValue();
    char * result = XMLString::transcode(value);
    return result;		
  }

  //Parse the attribute name and return its double value
  inline static double parseDouble(const char * name, DOMNamedNodeMap * map){						
    char * value = XMLString::transcode(getNode(name,map)->getNodeValue());
    double result = (double)atof(value);
    XMLString::release(&value);
    return result;
  }

  //returns the int value of the node
  inline static int parseInt(DOMNode * node){
    return XMLString::parseInt(node->getNodeValue());
  }

  //There are many nodes, but only one TEXT_NODE, this method
  //returns the string of this node
  inline static const char * parseContent(DOMNodeList * list){
    for (int i = 0, length = list->getLength(); i < length ; i++){
      DOMNode * node = list->item(i);
      if (node->getNodeType() == DOMNode::TEXT_NODE){
        return XMLString::transcode(node->getNodeValue());
      }
    }
    return "";
  }

  //Parse the attribute name and return its boolean value
  inline static bool parseBool(const char * name, DOMNamedNodeMap * map){
    DOMNode * node = getNode(name,map);
    char * nodeValue = XMLString::transcode(node->getNodeValue());
    bool result = XMLString::equals("true",nodeValue);		
    XMLString::release(&nodeValue);
    return result;
  }

  //Return if there is an attribute named name
  inline static bool contains(const char * name,DOMNamedNodeMap * map){
    XMLCh * translation = XMLString::transcode(name);
    bool result = map->getNamedItem(translation) != NULL;
    XMLString::release(&translation);
    return result;
  }

  //return if the element node has the same name as name
  inline static bool equals(DOMNode * node, const char * name){
    XMLCh * translation = XMLString::transcode(name);
    bool result = XMLString::equals(node->getNodeName(),translation);
    XMLString::release(&translation);
     return result;
  }
};

#endif
