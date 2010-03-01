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
    int ModelOptions::optSet;
    
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
		    options.set (PRINT_MODEL_OPTIONS);
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
		    cerr << "Unrecognised command-line option: --" << clo << endl << endl;
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
			options.set (PRINT_MODEL_OPTIONS);
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
	    << " -m --print-model	Print all model options with a non-default value and exit." << endl
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
    
    // Utility: converts option strings to codes and back
    class OptionCodeMap {
	// Lookup table to translate the strings used in the XML file to the internal enumerated values:
	map<string,OptionCodes> codeMap;
	
    public:
	OptionCodeMap () {
	    codeMap["PENALISATION_EPISODES"] = PENALISATION_EPISODES;
	    codeMap["NEGATIVE_BINOMIAL_MASS_ACTION"] = NEGATIVE_BINOMIAL_MASS_ACTION;
	    codeMap["ATTENUATION_ASEXUAL_DENSITY"] = ATTENUATION_ASEXUAL_DENSITY;
	    codeMap["LOGNORMAL_MASS_ACTION"] = LOGNORMAL_MASS_ACTION;
	    codeMap["NO_PRE_ERYTHROCYTIC"] = NO_PRE_ERYTHROCYTIC;
	    codeMap["MAX_DENS_CORRECTION"] = MAX_DENS_CORRECTION;
	    codeMap["INNATE_MAX_DENS"] = INNATE_MAX_DENS;
	    // 	codeMap["MAX_DENS_RESET"] = MAX_DENS_RESET;
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
	    codeMap["GARKI_DENSITY_BIAS"] = GARKI_DENSITY_BIAS;
	}
	
	OptionCodes operator[] (const string s) {
	    map<string,OptionCodes>::iterator codeIt = codeMap.find (s);
	    if (codeIt == codeMap.end()) {
		ostringstream msg;
		msg << "Unrecognised model option: ";
		msg << s;
		throw xml_scenario_error(msg.str());
	    }
	    return codeIt->second;
	}
	// reverse-lookup in map; only used for error/debug printing so efficiency is unimportant
	// doesn't ensure code is unique in the map either
	string toString (const OptionCodes code) {
	    for (map<string,OptionCodes>::iterator codeIt = codeMap.begin(); codeIt != codeMap.end(); ++codeIt) {
		if (codeIt->second == code)
		    return codeIt->first;
	    }
	    throw runtime_error ("toString called with unknown code");	// this is a code error
	}
    };
    
    void ModelOptions::init () {
	OptionCodeMap codeMap;
	
	// State of all default options:
	bitset<NUM_OPTIONS> defaultOptSet;
	defaultOptSet.set (MAX_DENS_CORRECTION);
	
	// Set optSet to defaults, then override any given in the XML file:
	optSet = defaultOptSet.to_ulong();
	
	const scnXml::OptionSet::OptionSequence& optSeq = InputData.getScenario().getModelOptions().getOption();
	for (scnXml::OptionSet::OptionConstIterator it = optSeq.begin(); it != optSeq.end(); ++it) {
	    if (it->getValue())
		optSet |= (1<<codeMap[it->getName()]);
	}
	
	cout << boolalpha;
	
	// Print non-default model options:
	if (CommandLine::option (CommandLine::PRINT_MODEL_OPTIONS)) {
	    cout << "Non-default model options:";
	    for (int i = 0; i < NUM_OPTIONS; ++i) {
		if (bool(optSet & (1<<i)) != defaultOptSet[i])
		    cout << "\t" << codeMap.toString(OptionCodes(i)) << "=" << bool(optSet & 1<<i);
	    }
	    cout << endl;
	}
	
	// Test for incompatible options
	
	// An incompatibility triangle, listing options incompatible with each option
	// Doesn't check required versions
	// Originally from "description of variables for interface" excel sheet
	bitset<NUM_OPTIONS> incompatibilities[NUM_OPTIONS];	// all default to 0
	
	incompatibilities[NEGATIVE_BINOMIAL_MASS_ACTION]	.set(LOGNORMAL_MASS_ACTION)
													.set(TRANS_HET)	.set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[LOGNORMAL_MASS_ACTION]		.set(TRANS_HET)	.set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)	.set(TRIPLE_HET);
	
	incompatibilities[ATTENUATION_ASEXUAL_DENSITY]	.set(INCLUDES_PK_PD)
													.set(DUMMY_WITHIN_HOST_MODEL)	.set(EMPIRICAL_WITHIN_HOST_MODEL);
	
	//Note: also MAX_DENS_CORRECTION is irrelevant when using DUMMY_WITHIN_HOST_MODEL or EMPIRICAL_WITHIN_HOST_MODEL
	// but we don't mark it incompatible so we can leave MAX_DENS_CORRECTION on by default.
	incompatibilities[INNATE_MAX_DENS]	.set(DUMMY_WITHIN_HOST_MODEL)	.set(EMPIRICAL_WITHIN_HOST_MODEL);
	
	incompatibilities[EMPIRICAL_WITHIN_HOST_MODEL]	.set(DUMMY_WITHIN_HOST_MODEL);
	
	incompatibilities[NON_MALARIA_FEVERS]	.set(MUELLER_PRESENTATION_MODEL);
	
	incompatibilities[TRANS_HET]			.set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)	.set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_HET]			.set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)	.set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[TREAT_HET]			.set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)	.set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_TRANS_HET]	.set(TRANS_TREAT_HET)	.set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[TRANS_TREAT_HET]		.set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_TREAT_HET]	.set(TRIPLE_HET);
	
	for (size_t i = 0; i < NUM_OPTIONS; ++i) {
	    if ((optSet & 1<<i) && (bitset<NUM_OPTIONS>(optSet) & incompatibilities[i]).any()) {
		ostringstream msg;
		msg << "Incompatible model options: " << codeMap.toString(OptionCodes(i)) << "=" << bool(optSet & 1<<i)
			<< " is incompatible with flags:";
		bitset<NUM_OPTIONS> incompat = (bitset<NUM_OPTIONS>(optSet) & incompatibilities[i]);
		for (int j = 0; j < NUM_OPTIONS; ++j) {
		    if (incompat[j])
			msg << "\t" << codeMap.toString(OptionCodes(i)) << "=" << bool(optSet & 1<<i);
		}
		throw xml_scenario_error (msg.str());
	    }
	}
	
	// Required options (above table can't check these):
	// NOTE: in principle, this restriction isn't required
	if ((optSet & EMPIRICAL_WITHIN_HOST_MODEL) && !(optSet & INCLUDES_PK_PD))
	    throw xml_scenario_error ("EMPIRICAL_WITHIN_HOST_MODEL requires INCLUDES_PK_PD");
	
	// Stop OM now if "--print-model" was on command-line
	if (CommandLine::option (util::CommandLine::PRINT_MODEL_OPTIONS))
	    throw cmd_exit ("Printed model version");
    }
    
} }