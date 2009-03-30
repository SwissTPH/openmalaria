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

#include <iostream>
#include <fstream>

using namespace std;

#include "boincWrapper.h"

#include "inputData.h" //Include parser for the input
#include "GSLWrapper.h" //Include wrapper for GSL library

#include "simulation.h"


/** main() - initializes and shuts down BOINC and GSL, loads scenario XML and
 * runs simulation. */
int main(int argc, char* argv[]){
  cout << "Message to std::cout"<<endl;
  try {
    string scenario_name = "scenario.xml";
    if (argc == 2) {
      scenario_name = argv[1];
    } else if (argc > 1) {
      cout << "Usage: "<<argv[0]<<" [file.xml]\nIf no file name is supplied, the name scenario.xml is used."<<endl;
      return 1;
    }
    
    BoincWrapper::init();
    
    scenario_name = BoincWrapper::resolveFile (scenario_name.c_str());
    
    //Change it and read it with boinc
    if (!createDocument(scenario_name)) {
      cerr << "APP. createDocument failed." <<endl;
      BoincWrapper::finish(-1);
      exit(-1);
    }
    
    GSL_SETUP();
    {
      Simulation simulation;	// CTOR runs
      simulation.start();
    }	// DTOR runs
    GSL_TEARDOWN();
    
    cleanDocument();
    BoincWrapper::finish(0);	// Never returns
  } catch (...) {
    cerr << "Error occurred." << endl;
    BoincWrapper::finish(-1);
  }
}
