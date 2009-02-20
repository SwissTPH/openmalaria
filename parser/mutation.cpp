
#ifndef MUTATION_CPP
#define MUTATION_CPP
#include "Names.h"

#include "converter.cpp"
#include "Constant.h"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Mutation : public MalariaNode {
private:
	int position;
public:
	Mutation(DOMNode * node){
		createNode(this,node);
	};
	
	~Mutation(void) {};
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
		position = MISSING_VALUE;
		bool hasPosition = Converter::contains(s_POSITION,map);
		if (hasPosition)
			position = Converter::parseInt(s_POSITION,map);		
	};

	void addChild(DOMNode * child){
		//Nothing to do
	}

	int getPosition(){		
		return position;
	}


#ifdef _LOG
	void debug(){
		cerr << "<mutation " ;
		cerr << "\tposition " << position <<
			"\t>\n";
	}
#endif

};

#endif
