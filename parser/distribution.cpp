/**
* Tiago Antao
*/

#ifndef DISTRIBUTION_CPP
#define DISTRIBUTION_CPP
/*test*/
#include "Names.h"

#include "proteomeInstance.cpp"

#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Distribution : public MalariaNode{
private:

	//The different proteomeInstances (called genome content to avoid name clashes)
	ProteomeInstance ** proteomeInstances;
	//The number of proteomeInstances
	int numProteomeInstances;
	int x;

public:

	Distribution(DOMNode * node):x(0){
		createNode(this,node);
	}

	~Distribution(void) {
		//We have to delete each gene
		for (int i = 0; i < numProteomeInstances; i++)
			delete proteomeInstances[i];
		//We can now delete the array
		delete [] proteomeInstances;
	}

	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		//Check the number of proteomeInstances
		numProteomeInstances = 0;
		for (int i =0, length = list->getLength(); i < length ; i++){
			DOMNode * child = list->item(i);
			//We have to check if it is really an element
			if (child->getNodeType() == DOMNode::ELEMENT_NODE)
				numProteomeInstances++;
		}
		//We can create the array
		proteomeInstances = new ProteomeInstance*[numProteomeInstances];
	};

	//Adds a new proteome instance: the method is called by its superclass
	//There can be only one element genome instances.
	void addChild(DOMNode * child){				
		proteomeInstances[x++] = new ProteomeInstance(child);
	}

	//return the number of proteome instances
	int getNumProteomeInstances(){
		return numProteomeInstances;
	}

	//return the (index)th proteome instance
	ProteomeInstance * getProteomeInstance(int index){
		return proteomeInstances[index];
	}

#ifdef _LOG
	void debug(){
		cerr << "<Distribution " <<	
			"\tnumProteomeInstances " << numProteomeInstances <<  
			"\t>\n";
	}
#endif
};

#endif
