/**
* Tiago Antao
*/

#ifndef PROTEOME_CPP
#define PROTEOME_CPP
/*test*/
#include "Names.h"

#include "content.cpp"
#include "distribution.cpp"

#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Proteome : public MalariaNode{
private:

	//Genome content
	Content * content;
	//Genome distribution
	Distribution * distribution;

public:

	Proteome(DOMNode * node){
		content = NULL;
		distribution = NULL;
		createNode(this,node);
	}

	~Proteome(void) {
		//We have to delete content and distribution

		delete content;
		delete distribution;
	}

	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
	}


	void addChild(DOMNode * child){				
		if (Converter::equals(child,s_CONTENT)){
			content = new Content(child);		
			return;
		}
		if (Converter::equals(child,s_DISTRIBUTION)){
			distribution = new Distribution(child);		
			return;
		}

	}

    Content* getContent() {
        return content;
    }

    Distribution* getDistribution() {
        return distribution;
    }

#ifdef _LOG
	void debug(){
		cerr << "<Proteome " <<  
			"\t>\n";
	}
#endif
};

#endif
