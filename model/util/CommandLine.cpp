/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "Global.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "util/BoincWrapper.h"
#include "util/StreamValidator.h"
#include "inputData.h"

#include <sstream>
#include <iostream>
#include <cassert>
#include <boost/lexical_cast.hpp>

namespace OM { namespace util {
    using boost::lexical_cast;
    
    bitset<CommandLine::NUM_OPTIONS> CommandLine::options;
    string CommandLine::resourcePath;
    double CommandLine::newEIR;
    string CommandLine::outputName;
    set<int> CommandLine::checkpoint_times;
    
    string parseNextArg (int argc, char* argv[], int& i) {
	++i;
	if (i >= argc)
	    throw runtime_error ("Expected an argument following the last option");
	return string(argv[i]);
    }
    
    string CommandLine::parse (int argc, char* argv[]) {
	options[COMPRESS_CHECKPOINTS] = true;	// turn on by default
	
	bool cloHelp = false, cloVersion = false, cloError = false;
	newEIR = numeric_limits<double>::quiet_NaN();
	bool fileGiven = false;
	string scenarioFile = "scenario.xml";
	outputName = "";
#	ifdef OM_STREAM_VALIDATOR
	string sVFile;
#	endif
	
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
		    if (fileGiven){
			throw runtime_error ("--scenario argument may only be given once");
		    }
		    scenarioFile = parseNextArg (argc, argv, i);
		    fileGiven = true;
		} else if (clo == "output") {
		    if (outputName != ""){
			throw runtime_error ("--output argument may only be given once");
		    }
		    outputName = parseNextArg (argc, argv, i);
		} else if (clo == "print-model") {
		    options.set (PRINT_MODEL_OPTIONS);
                    options.set (SKIP_SIMULATION);
		} else if (clo == "print-EIR") {
		    options.set (PRINT_ANNUAL_EIR);
                    options.set (SKIP_SIMULATION);
		} else if (clo == "set-EIR") {
		    if (newEIR == newEIR)	// initialised to NaN
			throw runtime_error ("--set-EIR already given");
		    newEIR = lexical_cast<double>(parseNextArg (argc, argv, i));
		    options.set (SET_ANNUAL_EIR);
                    options.set (SKIP_SIMULATION);
                } else if (clo == "validate-only") {
                    options.set (SKIP_SIMULATION);
		} else if (clo == "checkpoint") {
		    options.set (TEST_CHECKPOINTING);
		} else if (clo.compare (0,11,"checkpoint=") == 0) {
		    stringstream t;
		    t << clo.substr (11);
		    int time;
		    t >> time;
		    if (t.fail() || time <= 0) {
			cerr << "Expected: --checkpoint=t  where t is a positive integer" << endl;
			cloError = true;
			break;
		    }
		    checkpoint_times.insert( time );
		} else if (clo.compare (0,21,"compress-checkpoints=") == 0) {
		    stringstream t;
		    t << clo.substr (21);
		    bool b;
		    t >> b;
		    if (t.fail()) {
			cerr << "Expected: --compress-checkpoints=x  where x is 1 or 0" << endl;
			cloError = true;
			break;
		    }
		    options[COMPRESS_CHECKPOINTS] = b;
		} else if (clo == "checkpoint-duplicates") {
		    options.set (TEST_DUPLICATE_CHECKPOINTS);
#	ifdef OM_STREAM_VALIDATOR
		} else if (clo == "stream-validator") {
		    if (sVFile.size())
			throw runtime_error ("--stream-validator may only be given once");
		    sVFile = parseNextArg (argc, argv, i);
#	endif
                } else if (clo == "version") {
                    cloVersion = true;
                } else if (clo == "help") {
                    cloHelp = true;
		} else {
		    cerr << "Unrecognised command-line option: --" << clo << endl;
		    cloError = true;
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
			options.set (SKIP_SIMULATION);
		    } else if (clo[j] == 's') {
			if (fileGiven){
			    throw runtime_error ("-s argument may only be given once");
			}
			scenarioFile = parseNextArg (argc, argv, i);
			fileGiven = true;
		    } else if (clo[j] == 'o') {
			if (outputName != ""){
			    throw runtime_error ("-o argument may only be given once");
			}
			outputName = parseNextArg (argc, argv, i);
		    } else if (clo[j] == 'c') {
			options.set (TEST_CHECKPOINTING);
		    } else if (clo[j] == 'd') {
			options.set (TEST_DUPLICATE_CHECKPOINTS);
                    } else if (clo[j] == 'v') {
                        cloVersion = true;
                    } else if (clo[j] == 'h') {
                        cloHelp = true;
		    } else {
		    cerr << "Unrecognised command-line option: -" << clo[j] << endl;
		    cloError = true;
		    }
		}
	    } else {
		cerr << "Unexpected parameter: " << clo << endl << endl;
		cloError = true;
	    }
	}
	
	if (cloVersion || cloHelp){
            cerr<<"OpenMalaria simulator of malaria epidemiology and control, schema version "
                  <<InputDataType::SCHEMA_VERSION<<endl
                  <<"(oldest compatible: "<<InputDataType::SCHEMA_VERSION_OLDEST_COMPATIBLE
                  <<"). For more information, see"<<endl
                  <<"http://code.google.com/p/openmalaria/"<<endl<<endl
                  <<"OpenMalaria is copyright © 2005-2010 Swiss Tropical Institute and Liverpool"<<endl
                  <<"School Of Tropical Medicine."<<endl
                  <<"OpenMalaria comes with ABSOLUTELY NO WARRANTY. This is free software, and you"<<endl
                  <<"are welcome to redistribute it under certain conditions. See the file COPYING"<<endl
                  <<"or http://www.gnu.org/licenses/gpl-2.0.html for details of warranty or terms of"<<endl
                  <<"redistribution."<<endl<<endl;
        }
	if (cloHelp || cloError) {
	    cerr << "Usage: " << argv[0] << " [options]" << endl << endl
	    << "Options:"<<endl
	    << " -p --resource-path	Path to look up input resources with relative URLs (defaults to"<<endl
	    << "			working directory). Not used for output files."<<endl
	    << " -s --scenario file.xml	Uses file.xml as the scenario. If not given, scenario.xml is used." << endl
	    << "			If path is relative (doesn't start '/'), --resource-path is used."<<endl
	    << " -o --output file.txt	Uses file.txt as output file name. If not given, output.txt is used." << endl
	    << " -m --print-model	Print all model options with a non-default value and exit." << endl
	    << "    --print-EIR	Print the annual EIR (of each species in vector mode) and exit." << endl
	    << "    --set-EIR LEVEL	Scale the input EIR to a new annual level (innocs./person/year)"<<endl
	    << "			Note: updated XML file will be generated in working directory,"<<endl
	    << "			and will have other, mostly insignificant, differences to original."<<endl
            << "    --validate-only	Initialise and validate scenario, but don't run simulation." << endl
	    << "    --checkpoint=t	Forces a checkpoint a simulation time t. May be specified"<<endl
	    << "			more than once. Overrides --checkpoint option."<<endl
	    << " -c --checkpoint	Forces a checkpoint during each simulation"<<endl
	    << "			period, exiting after completing each"<<endl
	    << "			checkpoint. Doesn't require BOINC to do the checkpointing." <<endl
	    << " -d --checkpoint-duplicates"<<endl
	    << "			Write a checkpoint immediately after reading, which should be" <<endl
	    << "			identical to that read." <<endl
	    << "    --compress-checkpoints=boolean" << endl
	    << "			Set checkpoint compression on or off. Default is on." <<endl
#	ifdef OM_STREAM_VALIDATOR
	    << "    --stream-validator PATH" <<endl
	    << "			Use StreamValidator to validate against reference file PATH." <<endl
	    << "			(note: PATH must be absolute or relative to resource path)." <<endl
#	endif
            << " -v --version           Display the current schema version of OpenMalaria." << endl
            << " -h --help              Print this message." << endl<<endl
	    ;
	    if( cloError )
		throw std::invalid_argument( "bad argument" );
	}
	if( cloVersion || cloHelp ){
            throw cmd_exit ("Printed help");
        }
	
#	ifdef OM_STREAM_VALIDATOR
	if( sVFile.size() )
	    StreamValidator.loadStream( sVFile );
#	endif
	
	if (checkpoint_times.size())	// timed checkpointing overrides this
	    options[TEST_CHECKPOINTING] = false;

	if (outputName == ""){
	    outputName = "output.txt";
	}
#ifndef WITHOUT_BOINC
	outputName.append(".gz");
#endif

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
    
    /* These check parameters are as expected. They only really serve to make
     * sure important command-line parameters didn't change (and only in DEBUG mode)! */
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
    
} }
