/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "util/CommandLine.hpp"
#include "util/ModelOptions.hpp"
#include "util/errors.hpp"
#include "util/BoincWrapper.h"
#include "inputData.h"

#include <sstream>
#include <iostream>
#include <cassert>

/// Source code for CommandLine.hpp and ModelOptions.hpp
namespace OM { namespace util {
    
    bitset<CommandLine::NUM_OPTIONS> CommandLine::options;
    string CommandLine::resourcePath;
    std::bitset<NUM_OPTIONS> ModelOptions::optSet;
    
    string parseNextArg (int argc, char* argv[], int& i) {
	++i;
	if (i >= argc)
	    throw runtime_error ("Expected an argument following the last option");
	return string(argv[i]);
    }
    
    string CommandLine::parse (int argc, char* argv[]) {
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
		    if (resourcePath.size())
			throw runtime_error ("--resource-path (or -p) may only be given once");
		    resourcePath = parseNextArg (argc, argv, i).append ("/");
		} else if (clo == "scenario") {
		    if (fileGiven)
			throw runtime_error ("--scenario argument may only be given once");
		    scenarioFile = parseNextArg (argc, argv, i);
		    fileGiven = true;
		} else if (clo == "print-model") {
		    options.set (PRINT_MODEL_VERSION);
		} else if (clo == "checkpoint") {
		    options.set (TEST_CHECKPOINTING);
		} else if (clo.compare (0,21,"compress-checkpoints=") == 0) {
		    stringstream t;
		    t << clo.substr (21);
		    bool b;
		    t >> b;
		    options[COMPRESS_CHECKPOINTS] = b;
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
			if (resourcePath.size())
			    throw runtime_error ("--resource-path (or -p) may only be given once");
			resourcePath = parseNextArg (argc, argv, i).append ("/");
		    } else if (clo[j] == 'm') {
			options.set (PRINT_MODEL_VERSION);
		    } else if (clo[j] == 'c') {
			options.set (TEST_CHECKPOINTING);
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
    
    string CommandLine::lookupResource (const string& path) {
	string ret;
	if (path.size() >= 1 && path[0] == '/') {
	    // UNIX absolute path
	} else if (path.size() >= 3 && path[1] == ':'
	    && (path[2] == '\\' || path[2] == '/')) {
	    // Windows absolute path.. at least probably
	} else {	// relative
	    ret = resourcePath;
	}
	ret.append (path);
	return BoincWrapper::resolveFile (ret);
    }
    
    void CommandLine::staticCheckpoint (istream& stream) {
	string tOpt;
	string tResPath;
	tOpt & stream;
	tResPath & stream;
	assert (tOpt == options.to_string());
	assert (tResPath == resourcePath);
    }
    void CommandLine::staticCheckpoint (ostream& stream) {
	options.to_string() & stream;
	resourcePath & stream;
    }
    
    void ModelOptions::init () {
	// Lookup table to translate the strings used in the XML file to the internal enumerated values:
	map<string,OptionCodes> codeMap;
	codeMap["PENALISATION_EPISODES"] = PENALISATION_EPISODES;
	codeMap["NEGATIVE_BINOMIAL_MASS_ACTION"] = NEGATIVE_BINOMIAL_MASS_ACTION;
	codeMap["ATTENUATION_ASEXUAL_DENSITY"] = ATTENUATION_ASEXUAL_DENSITY;
	codeMap["LOGNORMAL_MASS_ACTION"] = LOGNORMAL_MASS_ACTION;
	codeMap["NO_PRE_ERYTHROCYTIC"] = NO_PRE_ERYTHROCYTIC;
	codeMap["MAX_DENS_CORRECTION"] = MAX_DENS_CORRECTION;
	codeMap["INNATE_MAX_DENS"] = INNATE_MAX_DENS;
	codeMap["MAX_DENS_RESET"] = MAX_DENS_RESET;
	codeMap["DUMMY_WITHIN_HOST_MODEL"] = DUMMY_WITHIN_HOST_MODEL;
	codeMap["PREDETERMINED_EPISODES"] = PREDETERMINED_EPISODES;
	codeMap["NON_MALARIA_FEVERS"] = NON_MALARIA_FEVERS;
	codeMap["INCLUDES_PK_PD"] = INCLUDES_PK_PD;
	codeMap["CLINICAL_EVENT_SCHEDULER"] = CLINICAL_EVENT_SCHEDULER;
	codeMap["MUELLER_PRESENTATION_MODEL"] = MUELLER_PRESENTATION_MODEL;
	codeMap["TRANS_HET"] = TRANS_HET;
	codeMap["COMORB_HET"] = COMORB_HET;
	codeMap["TREAT_HET"] = TREAT_HET;
	codeMap["COMORB_TRANS_HET"] = COMORB_TRANS_HET;
	codeMap["TRANS_TREAT_HET"] = TRANS_TREAT_HET;
	codeMap["COMORB_TREAT_HET"] = COMORB_TREAT_HET;
	codeMap["TRIPLE_HET"] = TRIPLE_HET;
	codeMap["EMPIRICAL_WITHIN_HOST_MODEL"] = EMPIRICAL_WITHIN_HOST_MODEL;
	
	// State of all default options:
	bitset<NUM_OPTIONS> defaultOptSet;
	defaultOptSet.set (MAX_DENS_CORRECTION);
	
	// Set optSet to defaults, then override any given in the XML file:
	optSet = defaultOptSet;
	
	const scnXml::ModelOptions::OptionSequence& optSeq = InputData.getScenario().getModelOptions().getOption();
	for (scnXml::ModelOptions::OptionConstIterator it = optSeq.begin(); it != optSeq.end(); ++it) {
	    map<string,OptionCodes>::iterator codeIt = codeMap.find (it->getName());
	    if (codeIt == codeMap.end()) {
		ostringstream msg;
		msg << "Unrecognised option: ";
		msg << it->getName();
		throw xml_scenario_error(msg.str());
	    }
	    optSet[codeIt->second] = it->getValue();
	}
	
	// Print non-default model options:
	if (CommandLine::option (CommandLine::PRINT_MODEL_VERSION)) {
	    cout << "Non-default model options:";
	    for (int i = 0; i < NUM_OPTIONS; ++i) {
		if (optSet[i] != defaultOptSet[i]) {
		    // reverse-lookup in map; only used for "--print-model" option so efficiency is unimportant:
		    for (map<string,OptionCodes>::iterator codeIt = codeMap.begin(); codeIt != codeMap.end(); ++codeIt) {
			if (codeIt->second == i) {
			    cout << "\t" << codeIt->first << "=" << optSet[i];
			    break;
			}
		    }
		}
	    }
	    cout << endl;
	}
	
	/* TODO: test incompatible options
	// Or'd flags of incompatibility triangle from
	// "description of variables for interface" excel sheet
	const int INCOMPATIBLITITIES[NUM_OPTIONS] = {
	    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD,	// PENALISATION_EPISODES
	    LOGNORMAL_MASS_ACTION | ANY_TRANS_HET,	// NEGATIVE_BINOMIAL_MASS_ACTION
	    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD | ANY_TRANS_HET,
	    ANY_TRANS_HET,
	    ANY_TRANS_HET,
	    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD | ANY_TRANS_HET,	// MAX_DENS_CORRECTION
	    DUMMY_WITHIN_HOST_MODEL | INCLUDES_PK_PD,
	    DUMMY_WITHIN_HOST_MODEL,
	    0,
	    0,
	    MUELLER_PRESENTATION_MODEL,
	    0,
	    0,
	    0,
	    COMORB_TRANS_HET | TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,
	    COMORB_TRANS_HET | TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,
	    COMORB_TRANS_HET | TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,
	    TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,
	    COMORB_TREAT_HET | TRIPLE_HET,
	    TRIPLE_HET,
	    0,
	    DUMMY_WITHIN_HOST_MODEL | (!INCLUDES_PK_PD)
	};
	
	for (size_t i = 0; i < NUM_VERSIONS; ++i) {
	    if (((optSet >> i) & 1) &&
		optSet & INCOMPATIBLITITIES[i]) {
		ostringstream msg;
	    msg << "Incompatible model versions: flag 1<<" << i
	    << " is incompatible with flags: ";
	    int incompat = (optSet & INCOMPATIBLITITIES[i]);
	    for (int i = 0; i < NUM_VERSIONS; ++i) {
		if ((incompat >> i) & 1)
		    msg << " 1<<"<<i;
	    }
	    //Note: this can occur if a version is listed as "incompatible with itself" in the above table
	    throw xml_scenario_error (msg.str());
	    }
	}
	*/
	if ((optSet[MAX_DENS_RESET]) && !(optSet[MAX_DENS_CORRECTION]))
	    throw xml_scenario_error ("MAX_DENS_RESET without MAX_DENS_CORRECTION doesn't make sense.");
	
	if (CommandLine::option (util::CommandLine::PRINT_MODEL_VERSION))
	    throw cmd_exit ("Printed model version");
    }
    
} }