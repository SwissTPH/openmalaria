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

#include <cstdlib>
#include "global.h"
#include "inputData.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

/*
Contains global variables and constants and utility functions that are used in different modules.

Constants (parameter)

*/

ModelVersion Global::modelVersion;
int Global::interval;
size_t Global::intervalsPerYear;
int Global::maxAgeIntervals;
int Global::simulationMode;
int Global::latentp;

CLO::CLO Global::clOptions = CLO::NONE;
bool Global::compressCheckpoints = true;

vector<int> Global::infantIntervalsAtRisk;
vector<int> Global::infantDeaths;

string Global::parseCommandLine (int argc, char* argv[]) {
  bool cloHelp = false;
  bool fileGiven = false;
  string scenarioFile = "scenario.xml";
  
  for (int i = 1; i < argc; ++i) {
    string clo = argv[i];
    
    // starts "--"
    if (clo.size() >= 2 && *clo.data() == '-' && clo.data()[1] == '-') {
      clo.assign (clo, 2, clo.size()-2);
      
      if (clo == "scenario") {
	if (!fileGiven) {
	  ++i;
	  if (i >= argc) {
	    cerr << "Expected: --scenario filename" << endl << endl;
	    cloHelp = true;
	    break;
	  }
	  scenarioFile = argv[i];
	  fileGiven = true;
	} else {
	  cerr << "Only one scenario file may be given." << endl;
	}
      } else if (clo == "print-model") {
	clOptions = CLO::CLO (clOptions | CLO::PRINT_MODEL_VERSION);
      } else if (clo == "enableERC") {
	clOptions = CLO::CLO (clOptions | CLO::ENABLE_ERC);
      } else if (clo == "checkpoint") {
	clOptions = CLO::CLO (clOptions | CLO::TEST_CHECKPOINTING);
      } else if (clo.compare (0,21,"compress-checkpoints=") == 0) {
	stringstream t;
	t << clo.substr (21);
	t >> compressCheckpoints;
	if (t.fail()) {
	  cerr << "Expected: --compress-checkpoints=x  where x is 1 or 0" << endl;
	  cloHelp = true;
	  break;
	}
      } else if (clo == "help") {
	cloHelp = true;
      } else {
	cerr << "Unrecognised option: --" << clo << endl << endl;
	cloHelp = true;
      }
    } else {
      cerr << "Unexpected parameter: " << clo << endl << endl;
      cloHelp = true;
    }
  }
  
  if (cloHelp) {
    cout << "Usage:" << endl
    << argv[0] << " [options]" << endl << endl
    << "Options:"<<endl
    << "  --scenario file.xml	Uses file.xml as the scenario. If not given, scenario.xml is used." << endl
    << "  --print-model		Print which model version flags are set (as binary with right-most"<<endl
    << "			digit representing option 0) and exit." << endl
    << "  --enableERC		Allow Emergence Rate Calculations (otherwise will stop if the"<<endl
    << "			values read from file are inaccurate)."<<endl
    << "  --checkpoint		Forces checkpointing once during simulation (during initialisation"<<endl
    << "			period), exiting after completing each"<<endl
    << "			checkpoint. Doesn't require BOINC to do the checkpointing." <<endl
    << "  --compress-checkpoints=boolean" << endl
    << "			Set checkpoint compression on or off. Default is on." <<endl
    << "  --help		Print this message." << endl
    ;
    exit (1);
  }
  
  return scenarioFile;
}

bool Global::initGlobal () {
  setModelVersion();
  interval=get_interval();
  if (daysInYear % interval !=  0) {
    cerr << "daysInYear not a multiple of interval" << endl;
    exit(-12);
  }
  intervalsPerYear = daysInYear/interval;
  infantDeaths.resize(intervalsPerYear);
  infantIntervalsAtRisk.resize(intervalsPerYear);
  latentp=get_latentp();
  maxAgeIntervals=(int)get_maximum_ageyrs()*intervalsPerYear;
  
  return clOptions & CLO::EARLY_EXIT;
}

void Global::setModelVersion () {
  modelVersion = (ModelVersion) get_model_version();
  
  if (clOptions & CLO::PRINT_MODEL_VERSION) {
    cout << "Model version: ";
    for (int i = NUM_VERSIONS; i >= 0; --i)
      cout << ((modelVersion >> i) & 1);
    cout << endl; 
  }
  
  // Or'd flags of incompatibility triangle from
  // "description of variables for interface" excel sheet
  const int INCOMPATIBLITITIES[NUM_VERSIONS] = {
    0,
    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD,	// 1
    LOGNORMAL_MASS_ACTION | LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM,	// 2
    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD | TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET,	// 3
    LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM,	// 4
    TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET,	// 5
    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD | TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET,	// 6
    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD,	// 7
    DUMMY_WITHIN_HOST_MODEL,	// 8
    0,	// 9
    0,	// 10
    MUELLER_PRESENTATION_MODEL,	// 11
    0,	// 12
    0,	// 13
    0,	// 14
    COMORB_TRANS_HET | TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,	// 15
    COMORB_TRANS_HET | TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,	// 16
    COMORB_TRANS_HET | TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,	// 17
    TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,	// 18
    COMORB_TREAT_HET | TRIPLE_HET,	// 19
    TRIPLE_HET,	// 20
    0,	// 21
  };
  
  for (size_t i = 0; i < NUM_VERSIONS; ++i)
    if (((modelVersion >> i) & 1) &&
          modelVersion & INCOMPATIBLITITIES[i]) {
      ostringstream msg;
      msg << "Incompatible model versions: flag " << i << " is incompatible with other flags.";
      throw xml_scenario_error (msg.str());
    }
  if (modelVersion & (MAX_DENS_CORRECTION | INNATE_MAX_DENS | MAX_DENS_RESET))
    cerr << "Warning: model version used is deprecated" << endl;
}

int Global::modIntervalsPerYear (int i) {
    int valmodIntervalsPerYear = i % intervalsPerYear;
    if (valmodIntervalsPerYear == 0) {
        valmodIntervalsPerYear=intervalsPerYear;
    }
    return valmodIntervalsPerYear;
}

xml_scenario_error::xml_scenario_error(const string&  __arg)
: runtime_error(__arg) { }
