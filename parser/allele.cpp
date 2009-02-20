
#ifndef ALLELE_CPP
#define ALLELE_CPP
#include "Names.h"

#include "converter.cpp"
#include "Constant.h"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Allele : public MalariaNode {
private:
    char* name;
	char* aminos;
	int CNV;
public:
	Allele(DOMNode * node){
		createNode(this,node);
	};
	
	~Allele(void) {};
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
        name = NULL;
		aminos = NULL;
        CNV = MISSING_VALUE;
		bool hasCNV = Converter::contains(s_CNV,map);
		if (hasCNV)
			CNV = Converter::parseInt(s_CNV,map);		
		bool hasAminos = Converter::contains(s_AMINOS,map);
		if (hasAminos)
			aminos = Converter::getValue(s_AMINOS,map);		
		bool hasName = Converter::contains(s_NAME,map);
		if (hasName)
			name = Converter::getValue(s_NAME,map);		
	};

	void addChild(DOMNode * child){
		//Nothing to do
	}

	int getCNV(){		
		return CNV;
	}

    char* getAminos() {
        return aminos;
    }

    char* getName() {
        return name;
    }


#ifdef _LOG
	void debug(){
		cerr << "<allele " ;
        if (name != NULL) {
		    cerr << "\tname " << name;
        }
        if (aminos != NULL) {
		    cerr << "\taminos " << aminos;
        }
        if (CNV != MISSING_VALUE) {
		    cerr << "\tCNV " << CNV;
        }
		cerr <<	"\t>\n";
	}
#endif

};

#endif
