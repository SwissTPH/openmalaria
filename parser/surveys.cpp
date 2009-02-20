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


#ifndef SURVEYS_CPP
#define SURVEYS_CPP

#include "Names.h"

#include "converter.cpp"
#include "intNode.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;


class Surveys : public MalariaNode
{
private:
	int * surveyTimes;	
	double detectionLimit;
	int summaryOption;
	int numGroups;
	int x;
	bool sort;
public:
	Surveys(DOMNode * node) :x(0){
		createNode(this,node);
	};

	~Surveys(void) {
		delete [] surveyTimes;
	};
	
	void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){	
		detectionLimit = Converter::parseDouble(s_DETECTION_LIMIT,map);
		summaryOption = Converter::parseInt(s_SUMMARY_OPTION,map);
		numGroups = 0;
		for (int i = 0, length = list->getLength(); i < length ; i++){
			DOMNode * child = list->item(i);
			if (child->getNodeType() == DOMNode::ELEMENT_NODE)
				numGroups++;
		}
		sort = true;
		surveyTimes = new int[numGroups];
	}

	void addChild(DOMNode * child){			
		IntNode time(child);
		surveyTimes[x++] = time.getValue();
		if (x > 1 && surveyTimes[x-2] > surveyTimes[x-1])
			sort = false;
	}

	double getDetectionLimit(){
		return detectionLimit;
	}

	int getSurvey(int index){
		return surveyTimes[index];
	}

	bool isSurvey(int time){
		if (sort){
			//Dichotomic search
			int minIndex = 0;
			int maxIndex = numGroups-1;
			if (surveyTimes[minIndex] == time) return true;
			if (surveyTimes[maxIndex] == time) return true;
			//We can decrement maxIndex and increment minIndex
			//Because we know that there are not at the border
			while (minIndex+1 < maxIndex){
				int middle = (maxIndex+minIndex)/2;
				int value = surveyTimes[middle];
				if (value == time) return true;
				if (value < time)
					maxIndex = middle;
				else 
					minIndex = middle;
			}
			return false;
		} else {
			//Brute force search
			for (int i = 0; i < numGroups; i++){
				if (time == surveyTimes[i])
					return true;
			}
			return false;
		}
	}

	int getSummaryOption(){
		return summaryOption;
	}

	int getNumSurveys(){
		return numGroups;
	}

#ifdef _LOG
	void debug(){
		cerr << "<surveys " << "\t" <<
			"\tsort " << sort <<
			"\tdetectionLimit " << detectionLimit << "\t" <<
			"summaryOption " << summaryOption << "\t" <<
			"numSurveys " << numGroups << "\t" <<
			">\n";
	}
#endif

};

#endif
