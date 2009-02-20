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


#ifndef MALARIANODE_H
#define MALARIANODE_H
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

/*
An abstract class for the nodes. 
*/
class MalariaNode
{	
public: 

  /**
   * A virtual destructor. To be C-ANSI compatible
   */
  virtual ~MalariaNode(void){};

	/**
	* Set the attributes and the list of node. Most of the time, it is useless to consider
	* the list: the method addChild will be called. 
	* For ageGroups, we need to know how many group there is. Then a loop is done on the list
	* to know the number of groups
	*
	*/
	virtual void setAttributes(DOMNamedNodeMap * map, DOMNodeList * list) = 0;

	/**
	* Only the nodes that are real nodes are called in this function
	*/
	virtual void addChild(DOMNode * child) = 0;

	//the debug function is used only if the 
#ifdef	_LOG
	/**
	* This function is called when we are in debug 
	mode and we want to see what was parsed
	*/
	virtual void debug() = 0;
#endif

protected:
	//
	// Create a new new node. It just simplify the use of the creation
	//
	inline void createNode(MalariaNode * node, DOMNode * domNode){
		//Get the list of internal entities
		DOMNodeList * list = domNode->getChildNodes();
		//set to the node the attributes and the list
		node->setAttributes(domNode->getAttributes(),list);

		//If we are in DEBUG mode we just print the node
#ifdef _LOG
		node->debug();
#endif

		DOMNode * child;		
		//Go through each node and give them to the object.
		for (int i = 0, length = list->getLength(); i < length; i++){
			child = list->item(i);
			if (child->getNodeType() == DOMNode::ELEMENT_NODE)
				node->addChild(child);
		}

	}
};

#endif
