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

/**@brief A wrapper around BOINC. This header should not include any BOINC headers!
 * 
 * The purpose of this is to avoid most other source files from having to
 * include BOINC headers, since these have caused some wierd issues with the
 * MSVC++ compiler. */

#ifndef Hmod_boincWrapper
#define Hmod_boincWrapper

#include <string>
using namespace std;

/// Wrapper around BOINC functions
namespace BoincWrapper {
  /// Initializes BOINC
  void init ();
  /// Cleans up BOINC and exits
  void finish (int err);
  
  /// Calls boinc_resolve_filename_s with inName, returning the result.
  string resolveFile (const string& inName);
  
  /// Report the proportion of work done
  void reportProgress (double progress);
  
  /// Returns true if it's time to write a checkpoint
  int timeToCheckpoint();
  /// Call when a checkpoint's completed
  void checkpointCompleted();
  
  /** Checks the current stream position (assumed to be end of file), resets
   * the stream, generates an MD5 sum from the stream (to EOF), and checks the
   * stream position is back where it started.
   * 
   * The idea of this is to generate a checksum of the file in a slightly
   * secure way by not closing and reopening the file. */
  void generateChecksum (istream& in);
}

/// Memory shared with graphics app:
namespace SharedGraphics {
  /// Creates shared memory and sets up updates.
  void init();
  
  const size_t KappaArraySize = 12;
  /// Function to set kappa in shared memory:
  void copyKappa(double *kappa);
}

#endif
