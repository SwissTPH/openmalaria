/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include <fstream>

using namespace std;

#include "BoincWrapper.h"

#include "inputData.h"	//Include parser for the input
#include "util/gsl.h"	//Include wrapper for GSL library

#include "global.h"
#include "simulation.h"


/** main() - initializes and shuts down BOINC and GSL, loads scenario XML and
 * runs simulation. */
int main(int argc, char* argv[]){
  int exitStatus = 0;
  try {
    string scenario_name =
      Global::parseCommandLine (argc, argv);
    
    BoincWrapper::init();
    
    scenario_name = BoincWrapper::resolveFile (scenario_name.c_str());
    
    //Change it and read it with boinc
    createDocument(scenario_name);
    
    Global::initGlobal();
    
    {
      Simulation simulation;	// constructor runs
      simulation.start();
    }	// simulation's destructor runs
  } catch (const cmd_exit& e) {	// this is not an error, but exiting due to command line
    cout << e.what() << "; exiting..." << endl;
  } catch (const exception& e) {
    cerr << "Exception: " << e.what() << endl;
    exitStatus = -1;
  } catch (...) {
    cerr << "Unknown exception" << endl;
    exitStatus = -1;
  }
  
  try {	// free XML memory (if allocated), and potentially save changes
    cleanDocument();
  } catch (...) {
    cerr << "cleanDocument failed" << endl;
    exitStatus = -1;
  }
  BoincWrapper::finish(exitStatus);	// Never returns
}
