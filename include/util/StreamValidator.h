/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef Hmod_StreamValidator
#define Hmod_StreamValidator

// Compile-time optional
#ifdef OM_STREAM_VALIDATOR
#include "Global.h"
#include <deque>
#include <boost/functional/hash.hpp>
#endif

namespace OM { namespace util {
    /// Use this function at validation points in code. If validator is not compile-
    /// time enabled, it will have no effect and should be optimised out.
    template<class T>
    inline void streamValidate (T x){
# ifdef OM_STREAM_VALIDATOR
	StreamValidator( x );
# endif
    }


// Compile-time optional
#ifdef OM_STREAM_VALIDATOR

    /** @brief StreamValidator implementation.
     *
     * Aim: help tracking down difficult checkpointing problems.
     * 
     * Usage:
     * 
     * 1. Enable cmake option OM_STREAM_VALIDATOR and compile.
     * 
     * 2. Run scenario as normal without checkpointing to generate a reference
     * trace (output to a file called "StreamValidator" in working directory).
     * 
     * 3. Run scenario with checkpointing and the "--stream-validator" option,
     * through a debugger. For example, from the build directory:
     * "test/run.py Empirical --gdb -- --stream-validator ../build/test/Empirical-6_Fsvp/StreamValidator --checkpoint"
     * In theory, the first section before checkpointing will be in-sync, so just
     * run ("run", "quit"). The second section after loading a checkpoint is
     * expected to run out of sync. Before running this, add a breakpoint at
     * the line mentioned in StreamValidator.cpp:
     * "break StreamValidator.cpp:103". Now run, and when the debugger stops
     * at your break-point, grab a stack trace ("bt"). If this doesn't accurately
     * enough show where the desync occurs, add some extra calls to the
     * "streamValidate()" macro in strategic code locations and repeat.
     */
    class StreamValidatorType {
    public:
	/// Create. Use store-mode unless loadStream() is called.
	StreamValidatorType () : storeMode(true) {}
	
	/// Save stream or confirm at end.
	void saveStream();
	
	/// Load a reference stream from file and switch to validation mode.
	void loadStream( const string& path );
	
	/// Templated function to take a value, hash it, and call handle.
	template<class T>
	void operator() (T value){
	    // creating a copy for each use: probably not efficient
	    // can provide specialized templates if necessary
	    boost::hash<T> hash_func;
	    
	    handle( hash_func( value ) );
	}
	
	/// Either store in reference stream or validate against reference stream.
	void handle( size_t value );
	
	/// Checkpointing
	template<class S>
	void operator& (S& cp_str) {
	    checkpoint( cp_str );
	}
	
    private:
	void checkpoint ( istream& cp_str );
	void checkpoint ( ostream& cp_str ) const;
	
	// True: read and store a reference value stream.
	// False: validate against a reference stream.
	bool storeMode;
	
	// Next position in the stream to read from (validation mode only).
	deque<size_t>::iterator readIt;
	
	// Essentially we want to store a massive sequence of data in-memory
	// (simplist with regard to checkpointing.) std::deque is probably an
	// efficient way of doing this.
	deque<size_t> stream;
    };
    extern StreamValidatorType StreamValidator;

#endif	//OM_STREAM_VALIDATOR
} }
#endif	//Hmod_StreamValidator
