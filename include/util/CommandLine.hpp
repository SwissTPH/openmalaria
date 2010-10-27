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
#include <set>
#include <bitset>
#include <limits>
using namespace std;

namespace OM { namespace util {
    /// Command line options and processing
    class CommandLine {
    public:
	/// Boolean command-line options
	enum Options {
	    /** Outputs non-default "ModelOptions" values in a human-readable form. */
	    PRINT_MODEL_OPTIONS = 0,
	    /** Forces checkpoints in the middle of each simulation phase,
	     * exiting immediately afterwards.
	     * 
	     * Also see checkpoint_times, which overrides this option. */
	    TEST_CHECKPOINTING,
	    /** Write a checkpoint immediately after loading one. Allows
	     * confirmation that a duplicate is produced. */
	    TEST_DUPLICATE_CHECKPOINTS,
	    /** Compress checkpoint files with gzip before writing.
	     * Even with binary checkpoints, this has a big effect. */
	    COMPRESS_CHECKPOINTS,
	    /** Do initialisation and error checks, but don't run simulation. */
	    SKIP_SIMULATION,
	    /** Print the annual EIR. */
	    PRINT_ANNUAL_EIR,
	    /** Scale EIR to a new annual level. */
	    SET_ANNUAL_EIR,
	    NUM_OPTIONS
	};
	
	/** Return true if given option (from CommandLine::Options) is active. */
	static inline bool option(size_t code) {
	    return options.test(code);
	}
	
	/** Return first checkpointing timestep _greater than_ timestep passed,
	 * or min int value if no (more) checkpoint times. */
	static int getNextCheckpointTime( int now ) {
	    set<int>::iterator it = checkpoint_times.upper_bound( now );
	    if( it == checkpoint_times.end() )
		return numeric_limits< int >::min();
	    else
		return *it;
	}
	
	/** Prepend if path is relative, prepend it with clResourcePath.
	* Then passes the resulting (or original) path through
	* BoincWrapper::resolveFile() and returns the result. */
	static string lookupResource (const string& path);
	
	/** Get the new target annual EIR (only set if SET_ANNUAL_EIR option
	 * is true). */
	static inline double getNewEIR (){
	    return newEIR;
	}
	
	/** Get the name of the output file.
	 */
	static inline string getOutputName (){
	    return outputName;
	}

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
	// Static parameters: don't require checkpointing (set by parse())
	static bitset<NUM_OPTIONS> options;	// boolean options
	static string resourcePath;
	
	static double newEIR;
	//Output filename (for main output file "output.txt")
	static string outputName;

	/** Set of simulation times at which a checkpoint should be written and
	* program should exit (to allow resume). */
	static set<int> checkpoint_times;
    };
} }
#endif
