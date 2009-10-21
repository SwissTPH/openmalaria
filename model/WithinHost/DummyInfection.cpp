/*
 This file is part of OpenMalaria.

 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

*/

#include "WithinHost/DummyInfection.h"
#include "inputData.h"
#include "util/gsl.h"
#include <algorithm>
#include <sstream>
#include <string.h>


DummyInfection::~DummyInfection() {
}
void DummyInfection::destroy() {
}

void DummyInfection::init (){
}

DummyInfection::DummyInfection(int simulationTime) :
  Infection (simulationTime)
{
    //Initialize current infection data
    _startdate=simulationTime;
    _density=8;	// increased by DH to avoid zeros in initialKappa
    _duration=infectionDuration(); 
    
    if (Global::modelVersion & INCLUDES_PK_PD)
      _proteome = ProteomeInstance::newInfection();
    else
      _proteome = NULL;
}

int DummyInfection::getEndDate(){
  return _startdate+_duration/Global::interval;
}

int DummyInfection::infectionDuration(){
    //We set some arbitrary maximum duration.
    return 100;
}

void DummyInfection::write (ostream& out) const {
  writeInfection (out);
  out << _duration << endl; 
}

DummyInfection::DummyInfection (istream& in) :
  Infection (in)
{
  in >> _duration;
}

void DummyInfection::determineWithinHostDensity(){
  const double GROWTH_RATE = 8.0;
  const double PARASITE_THRESHOLD = 1;
  
  /*
    if the density gets to be < 1 parasite per host then clear infections
    is called by making the duration negative. 
    */
  if (_density < PARASITE_THRESHOLD) {
    _duration=-99;
    _density = 0.0;
  } else {
    _density = (int(_density*GROWTH_RATE) % 20000);
  }
}
