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

#ifndef Hmod_util_CommandLine
#define Hmod_util_CommandLine

#include <string>
#include <bitset>
using namespace std;

namespace OM { namespace util {
    /// Command line options and processing
    class CommandLine {
    public:
	/// Boolean command-line options
	enum Options {
	    PRINT_MODEL_VERSION = 0,	// outputs modelVersion in a human-readable form
	    TEST_CHECKPOINTING,	// forces a checkpoint in the middle of initialisation, followed by exiting
	    COMPRESS_CHECKPOINTS,	// compress checkpoint files with gzip before writing
	    NUM_OPTIONS
	};
	
	/** Return true if given option (from CommandLine::Options) is active. */
	static inline bool option(size_t code) {
	    return options.test(code);
	}
	/** Prepend if path is relative, prepend it with clResourcePath.
	* Then passes the resulting (or original) path through
	* BoincWrapper::resolveFile() and returns the result. */
	static string lookupResource (const string& path);
	
	/** Looks through all command line options.
	*
	* @returns The name of the scenario XML file to use.
	*
	* Throws cmd_exit in the case a help message is printed. Help
	* is printed to cout, which necessitate calling this function
	* before BOINC is initialised.
	* 
	* In other cases command-line parameters cause variables to be set in Global
	* to achieve the desired result. */
	static string parse (int argc, char* argv[]);
	
	/** @brief Checkpointing.
	*
	* Not really required; more to confirm things are expected. */
	static void staticCheckpoint (istream& stream);
	static void staticCheckpoint (ostream& stream);	///< ditto
	
    private:
	// Static parameters: set by parse _and_ checked on checkpoint load
	static bitset<NUM_OPTIONS> options;	// boolean options
	static string resourcePath;
    };
} }
#endif
