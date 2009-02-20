/**
* Tiago Antao
*/

#ifndef GENE_CPP
#define GENE_CPP
/*test*/
#include "Names.h"

#include "mutation.cpp"

#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Protein : public MalariaNode{
private:

	//The protein name
	char * name;	
	//The different proteins
	Mutation ** mutations;
	//The number of proteins
	int numMutations;
	int x;

public:

	Protein(DOMNode * node):x(0){
		createNode(this,node);
	}

	~Protein(void) {
        delete name;
		//We have to delete each group
		for (int i = 0; i < numMutations; i++)
			delete mutations[i];
		//We can now delete the array
		delete [] mutations;
	}

	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		//Get the name
		name =  Converter::getValue(s_NAME,map);				
		//Check the number of mutations
		numMutations = 0;
		for (int i =0, length = list->getLength(); i < length ; i++){
			DOMNode * child = list->item(i);
			//We have to check if it is really an element
			if (child->getNodeType() == DOMNode::ELEMENT_NODE)
				numMutations++;
		}
		//We can create the array
		mutations = new Mutation*[numMutations];
	};

	//Adds a new group: the method is called by its superclass
	//There can be only one element group.
	void addChild(DOMNode * child){				
		mutations[x++] = new Mutation(child);
	}

	//Return the name
	char* getName(){
		return name;
	}

	//return the number of mutations
	int getNumMutations(){
		return numMutations;
	}

	//return the (index)th mutation
	Mutation * getMutation(int index){
		return mutations[index];
	}

#ifdef _LOG
	void debug(){
		cerr << "<Gene " <<  
			"\tname " << name << 
			"\tnumMutations " << numMutations <<  
			"\t>\n";
	}
#endif
};

#endif
