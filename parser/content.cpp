/**
* Tiago Antao
*/

#ifndef CONTENT_CPP
#define CONTENT_CPP
/*test*/
#include "Names.h"

#include "protein.cpp"

#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Content : public MalariaNode{
private:

	//The different genes 
	Protein ** proteins;
	//The number of Proteins
	int numProteins;
	int x;

public:

	Content(DOMNode * node):x(0){
		createNode(this,node);
	}

	~Content(void) {
		//We have to delete each gene
		for (int i = 0; i < numProteins; i++)
			delete proteins[i];
		//We can now delete the array
		delete [] proteins;
	}

	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		//Check the number of genes
		numProteins = 0;
		for (int i =0, length = list->getLength(); i < length ; i++){
			DOMNode * child = list->item(i);
			//We have to check if it is really an element
			if (child->getNodeType() == DOMNode::ELEMENT_NODE)
				numProteins++;
		}
		//We can create the array
		proteins = new Protein *[numProteins];
	};

	//Adds a new proteins: the method is called by its superclass
	//There can be only one element proteins.
	void addChild(DOMNode * child){				
		proteins[x++] = new Protein(child);
	}

	//return the number of Protein
	int getNumProteins(){
		return numProteins;
	}

	//return the (index)th Protein
	Protein * getProtein(int index){
		return proteins[index];
	}

#ifdef _LOG
	void debug(){
		cerr << "<Content " <<	
			"\tnumProteins " << numProteins <<  
			"\t>\n";
	}
#endif
};

#endif
