/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

/**@brief A wrapper around BOINC. This header should not include any BOINC headers!
 * 
 * The purpose of this is to avoid most other source files from having to
 * include BOINC headers, since these have caused some wierd issues with the
 * MSVC++ compiler. */

#ifndef Hmod_boincWrapper
#define Hmod_boincWrapper

#include <string>

namespace OM { namespace util {
/// Wrapper around BOINC functions
namespace BoincWrapper {
  /// Initializes BOINC
  void init ();
  /// Cleans up BOINC and exits
  void finish (int err);
  
  /** Calls boinc_resolve_filename_s with inName, returning the result.
   * Needs to be used for input and output files. */
  std::string resolveFile (const std::string& inName);
  
  /// Check whether a file exists. If using BOINC, use it's filesystem function;
  /// otherwise use an istream.
  bool fileExists(const char* path);
  
  /// Report the proportion of work done
  /// (progress is now / duration)
  void reportProgress (int now, int duration);
  
  /// Returns true if it's time to write a checkpoint
  int timeToCheckpoint();
  /// Call when a checkpoint's completed
  void checkpointCompleted();
  
  /// Open a critical section (see http://boinc.berkeley.edu/trac/wiki/BasicApi)
  void beginCriticalSection();
  /// End a critical section
  void endCriticalSection();
}

/** Class for storing and generating checksums.
*
* These are nearly md5 sums but not quite.
**********************************************************/
class Checksum {
public:
    /** Checks the current stream position (assumed to be end of file), resets
    * the stream, generates a checksum from the stream (to EOF), and checks the
    * stream position is back where it started.
    * 
    * Checksum generated is a (potentially perturbed) MD5 sum.
    * 
    * The idea of this is to generate a checksum of the file in a slightly
    * secure way by not closing and reopening the file. */
    static Checksum generate (std::istream& fileStream);
    
    /// Copy ctor to make sure it copies data, not pointers
    Checksum (const Checksum& that) {
	for(int i = 0; i < 16; ++i)
	    data[i] = that.data[i];
    }
    
    /** Write the checksum data to a file in base 16. */
    void writeToFile (std::string filename);
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
	for(int i = 0; i < 16; ++i)
	    data[i] & stream;
    }
    
    bool operator!= (Checksum& that) {
	for(int i = 0; i < 16; ++i)
	    if (data[i] != that.data[i])
		return true;
	return false;
    }
    
private:
    Checksum () {}		// use static functions to create
    unsigned char data[16];	// checksum
};

} }
#endif
