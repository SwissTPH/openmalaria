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

#ifndef ANOPHELES_CPP
#define ANOPHELES_CPP

#include "Names.h"

#include "converter.cpp"

#include "MalariaNode.h"
#include <iostream>
using namespace std;

class EIR : public MalariaNode{
public:
  double a0;
  double a1;
  double b1;
  double a2;
  double b2;
  double EIRRotateAngle;
 
  EIR(DOMNode * node){
    createNode(this,node);
  };
 
  void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
    a0 = Converter::parseDouble("a0",map);
    a1 = Converter::parseDouble("a1",map);
    a2 = Converter::parseDouble("a2",map);
    b1 = Converter::parseDouble("b1",map);
    b2 = Converter::parseDouble("b2",map);
    EIRRotateAngle = Converter::parseDouble("EIRRotateAngle",map);
  }
  void addChild(DOMNode * child){
  }
};

class Mosq : public MalariaNode{
public:
  int mosqRestDuration;
  double mosqSeekingDeathRate;
  double mosqSeekingDuration;
  double mosqProbBiting;
  double mosqProbFindRestSite;
  double mosqProbResting;
  double mosqProbOvipositing;

  Mosq(DOMNode * node){
    createNode(this,node);
  };

  void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
    mosqRestDuration = Converter::parseInt("mosqRestDuration",map);
    mosqSeekingDeathRate = Converter::parseDouble("mosqSeekingDeathRate",map);
    mosqSeekingDuration = Converter::parseDouble("mosqSeekingDuration",map);
    mosqProbBiting = Converter::parseDouble("mosqProbBiting",map);
    mosqProbFindRestSite = Converter::parseDouble("mosqProbFindRestSite",map);
    mosqProbResting = Converter::parseDouble("mosqProbResting",map);
    mosqProbOvipositing = Converter::parseDouble("mosqProbOvipositing",map);
  }
  void addChild(DOMNode * child){
  }

};

class Anopheles : public MalariaNode{

private:	
  EIR* eir;
  Mosq* mosq;
  bool useNv0Guess;

public:
  Anopheles(DOMNode * node){
    createNode(this,node);
  };

  ~Anopheles(void){
  };

  void setAttributes(DOMNamedNodeMap * map,DOMNodeList * list){
    useNv0Guess = Converter::parseBool(s_USENV0GUESS,map);
  }

  void addChild(DOMNode * child){	
    if (Converter::equals(child,s_EIR)){
      eir = new EIR(child);
    }
    if (Converter::equals(child,s_ANOPHELES)){
      mosq = new Mosq(child);
    }
  }

  bool isUseNv0Guess(){
    return useNv0Guess;
  }
  double geta0(){
    return eir->a0;
  }

  double geta1(){
    return eir->a1;
  }

  double geta2(){
    return eir->a2;
  }

  double getb1(){
    return eir->b1;
  }

  double getb2(){
    return eir->b2;
  }

  double getEIRRotateAngle(){
    return eir->EIRRotateAngle;
  }

  double getMosqRestDuration(){
    return mosq->mosqRestDuration;
  }

  double getMosqSeekingDeathRate(){
    return mosq->mosqSeekingDeathRate;
  }
  double getMosqSeekingDuration(){
    return mosq->mosqSeekingDuration;
  }
  double getMosqProbBiting(){
    return mosq->mosqProbBiting;
  }
  double getMosqProbFindRestSite(){
    return mosq->mosqProbFindRestSite;
  }

  double getMosqProbResting(){
    return mosq->mosqProbResting;
  }

  double getMosqProbOvipositing(){
    return mosq->mosqProbOvipositing;
  }

#ifdef _LOG
	void debug(){
	}
#endif

};

#endif
