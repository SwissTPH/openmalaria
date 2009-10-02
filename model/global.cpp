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
#include <cmath>
#include <stdexcept>

#include "global.h"
#include "inputData.h"
#include "util/BoincWrapper.h"

/*
Contains global variables and constants and utility functions that are used in different modules.

Constants (parameter)

*/

ModelVersion Global::modelVersion;
int Global::interval;
size_t Global::intervalsPerYear;
double Global::yearsPerInterval;
int Global::maxAgeIntervals;

CLO::CLO Global::clOptions = CLO::NONE;
string Global::clResourcePath;
bool Global::compressCheckpoints = true;

string parseNextArg (int argc, char* argv[], int& i) {
  ++i;
  if (i >= argc)
    throw runtime_error ("Expected an argument following the last option");
  return string(argv[i]);
}

string Global::parseCommandLine (int argc, char* argv[]) {
  bool cloHelp = false;
  bool fileGiven = false;
  string scenarioFile = "scenario.xml";
  
  /* Simple command line parser. Seems to work fine.
   * If an extension is wanted, http://tclap.sourceforge.net/ looks good. */
  for (int i = 1; i < argc; ++i) {
    string clo = argv[i];
    
    // starts "--"
    if (clo.size() >= 2 && *clo.data() == '-' && clo.data()[1] == '-') {
      clo.assign (clo, 2, clo.size()-2);
      
      if (clo == "resource-path") {
	if (clResourcePath.size())
	  throw runtime_error ("--resource-path (or -p) may only be given once");
	clResourcePath = parseNextArg (argc, argv, i).append ("/");
      } else if (clo == "scenario") {
	if (fileGiven)
	  throw runtime_error ("--scenario argument may only be given once");
	scenarioFile = parseNextArg (argc, argv, i);
	fileGiven = true;
      } else if (clo == "print-model") {
	clOptions = CLO::CLO (clOptions | CLO::PRINT_MODEL_VERSION);
      } else if (clo == "enableERC") {
	clOptions = CLO::CLO (clOptions | CLO::ENABLE_ERC);
      } else if (clo == "noErcValidation") {
	clOptions = CLO::CLO (clOptions | CLO::NO_ERC_VALIDATION);
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
    } else if (clo.size() >= 1 && *clo.data() == '-') {	// single - (not --)
      for (size_t j = 1; j < clo.size(); ++j) {
	if (clo[j] == 'p') {
	  if (j + 1 != clo.size())
	    throw runtime_error ("a path must be given as next argument after -p");
	  if (clResourcePath.size())
	    throw runtime_error ("--resource-path (or -p) may only be given once");
	  clResourcePath = parseNextArg (argc, argv, i).append ("/");
	} else if (clo[j] == 'v') {
	  clOptions = CLO::CLO (clOptions | CLO::NO_ERC_VALIDATION);
	} else if (clo[j] == 'm') {
	  clOptions = CLO::CLO (clOptions | CLO::PRINT_MODEL_VERSION);
	} else if (clo[j] == 'c') {
	  clOptions = CLO::CLO (clOptions | CLO::TEST_CHECKPOINTING);
	} else if (clo[j] == 'h') {
	  cloHelp = true;
	}
      }
    } else {
      cerr << "Unexpected parameter: " << clo << endl << endl;
      cloHelp = true;
    }
  }
  
  if (cloHelp) {
    cout << "Usage: " << argv[0] << " [options]" << endl << endl
    << "Options:"<<endl
    << " -p --resource-path	Path to look up input resources with relative URLs (defaults to"<<endl
    << "			working directory). Not used for output files."<<endl
    << "    --scenario file.xml	Uses file.xml as the scenario. If not given, scenario.xml is used." << endl
    << "			If path is relative (doesn't start '/'), --resource-path is used."<<endl
    << " -m --print-model	Print which model version flags are set (as binary with right-most"<<endl
    << "			digit representing option 0) and exit." << endl
    << "    --enableERC		Allow Emergence Rate Calculations (otherwise will stop if the"<<endl
    << "			values read from file are inaccurate)."<<endl
    << " -v --noErcValidation	Disables most of the emergence rate validation when human"<<endl
    << "			infectiousness validates."<<endl
    << " -c --checkpoint	Forces checkpointing once during simulation (during initialisation"<<endl
    << "			period), exiting after completing each"<<endl
    << "			checkpoint. Doesn't require BOINC to do the checkpointing." <<endl
    << "    --compress-checkpoints=boolean" << endl
    << "			Set checkpoint compression on or off. Default is on." <<endl
    << " -h --help		Print this message." << endl
    ;
    throw cmd_exit ("Printed help");
  }
  
  return scenarioFile;
}

void Global::initGlobal () {
  setModelVersion();
  interval=get_interval();
  if (daysInYear % interval !=  0) {
    cerr << "daysInYear not a multiple of interval" << endl;
    exit(-12);
  }
  intervalsPerYear = daysInYear/interval;
  yearsPerInterval = double(interval) / double(daysInYear);
  maxAgeIntervals=(int)get_maximum_ageyrs()*intervalsPerYear;
}

string Global::lookupResource (const string& path) {
  string ret;
  if (path.size() >= 1 && path[0] == '/') {
    // UNIX absolute path
  } else if (path.size() >= 3 && path[1] == ':'
	  && (path[2] == '\\' || path[2] == '/')) {
    // Windows absolute path.. at least probably
  } else {	// relative
    ret = clResourcePath;
  }
  ret.append (path);
  return BoincWrapper::resolveFile (ret);
}

void Global::setModelVersion () {
  modelVersion = (ModelVersion) get_model_version();
  
  // Or'd flags of incompatibility triangle from
  // "description of variables for interface" excel sheet
  const int INCOMPATIBLITITIES[NUM_VERSIONS] = {
    1,	// non-existent, so make it incompatible with itself
    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD,	// 1
    LOGNORMAL_MASS_ACTION | ANY_TRANS_HET,	// 2
    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD | ANY_TRANS_HET,	// 3
    ANY_TRANS_HET,	// 4
    ANY_TRANS_HET,	// 5
    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD | ANY_TRANS_HET,	// 6
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
    TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,   // 18
    COMORB_TREAT_HET | TRIPLE_HET,	// 19
    TRIPLE_HET,	// 20
    0,	// 21
    DUMMY_WITHIN_HOST_MODEL | (!INCLUDES_PK_PD) // 22
  };
  
  for (size_t i = 0; i < NUM_VERSIONS; ++i)
    if (((modelVersion >> i) & 1) &&
          modelVersion & INCOMPATIBLITITIES[i]) {
      ostringstream msg;
      msg << hex << "Incompatible model versions: flag 0x" << (1<<i)
	  << " is incompatible with flags: 0x" << (modelVersion & INCOMPATIBLITITIES[i]);
      //Note: this can occur if a version is listed as "incompatible with itself" in the above table
      throw xml_scenario_error (msg.str());
    }
  if (modelVersion & (MAX_DENS_CORRECTION | INNATE_MAX_DENS | MAX_DENS_RESET))
    cerr << "Warning: model version used is deprecated" << endl;
  
  if (clOptions & CLO::PRINT_MODEL_VERSION) {
    cout << "Model flags:";
    for (int i = 0; i < NUM_VERSIONS; ++i) {
      if ((modelVersion >> i) & 1)
	cout << " 1<<"<<i;
    }
    cout << endl;
    throw cmd_exit ("Printed model version");
  }
}

void Global::validateListSize (long length) {
  if (length < 0 || length > 1000) {
    ostringstream s;
    s << "List length out of range: " << length;
    throw checkpoint_error(s.str());
  }
}

xml_scenario_error::xml_scenario_error(const string&  __arg)
  : runtime_error(__arg) { }

checkpoint_error::checkpoint_error(const string&  __arg)
  : runtime_error(string("Error reading checkpoint: ").append(__arg)) { }

cmd_exit::cmd_exit(const string& __arg)
  : runtime_error(__arg) { }
