/**
* Tiago Antao
*/

#ifndef PROTEOME_INSTANCE_CPP
#define PROTEOME_INSTANCE_CPP
#include "Names.h"

#include "allele.cpp"

#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class ProteomeInstance : public MalariaNode{
private:
    //proportion (100.0 is max)
    double proportion;

    //fitness (100.0 would be wild-type normally)
	double fitness;

	//Alleles for each gene
	Allele ** alleles;
	//The number of alleles
	int numAlleles;
	int x;

public:

	ProteomeInstance(DOMNode * node):x(0){
		createNode(this,node);
	}

	~ProteomeInstance(void) {
		//We have to delete each allele
		for (int i = 0; i < numAlleles; i++)
			delete alleles[i];
		//We can now delete the array
		delete [] alleles;
	}

	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		//Get the proportion: the value is a double
        proportion =  Converter::parseDouble(s_PROPORTION,map);
		//Get the fitness: the value is a double
        fitness =  Converter::parseDouble(s_FITNESS,map);
		//Check the number of alleles
		numAlleles = 0;
		for (int i =0, length = list->getLength(); i < length ; i++){
			DOMNode * child = list->item(i);
			//We have to check if it is really an element
			if (child->getNodeType() == DOMNode::ELEMENT_NODE)
				numAlleles++;
		}
		//We can create the array
		alleles = new Allele *[numAlleles];
	};

	//Adds a new allele: the method is called by its superclass
	void addChild(DOMNode * child){				
		alleles[x++] = new Allele (child);
	}

	//Return the proportion
	double getProportion(){
		return proportion;
	}

	//Return the fitness
	double getFitness(){
		return fitness;
	}

	//return the number of alleles
	int getNumAlleles(){
		return numAlleles;
	}

	//return the (index)th alleles
	Allele * getAllele(int index){
		return alleles[index];
	}

#ifdef _LOG
	void debug(){
		cerr << "<ProteomeInstance " <<  
			"\tproportion " << proportion << 
			"\tfitness " << fitness << 
			"\tnumAlleles " << numAlleles <<  
			"\t>\n";
	}
#endif
};

#endif
