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
#ifndef Hmod_Global
#define Hmod_Global

#ifdef _MSC_VER
#pragma warning(disable: 4290)	// disable some warnings on MSVC
#define _CRT_SECURE_NO_DEPRECATE
#define finite(x) _finite(x)
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <fcntl.h>
#include <math.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>
using namespace std;

#include "Constant.h"
#include "util/errors.hpp"
#include "util/checkpoint.hpp"
using namespace OM::util::checkpoint;


/// Command-Line Option possibilities
namespace CLO {
  enum CLO {
    /* Values here are written in hexadecimal: http://en.wikipedia.org/wiki/Hexadecimal
     * Many are designed to be "flags", so the value corresponds to a single bit:
     * http://en.wikipedia.org/wiki/Flag_byte */
    NONE		= 0x0,
    
    PRINT_MODEL_VERSION	= 0x1,	// outputs modelVersion in a human-readable form
    TEST_CHECKPOINTING	= 0x2,	// forces a checkpoint in the middle of initialisation, followed by exiting
  };
}

class Global
{
public:
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
  static string parseCommandLine (int argc, char* argv[]);
  
  /** Sets parameters in Global and performs some checks.
   *
   * Also outputs some extra information for command-line options.
   * 
   * @throws cmd_exit if we should stop before running the simulation because
   * of an option. */
  static void initGlobal ();
  
  /// Variables that must be checkpointed.
  //@{

  /** Model version defines which implementations of hard-coded options should be
   * used. The integer value of modelVersion passed from the .xml is converted to
   * binary with each bit corresponding to a different dichotomous option.  The
   * original default model is modelVersion=0 */
  static ModelVersion modelVersion;
  //@}
  
  /// Data read from xml which doesn't need to be checkpointed.
  //@{

  /// temporal resolution of simulation, in days
  static int interval;
  /// Number of timesteps in 5 days
  static int intervalsPer5Days;
   //Simulation time steps per year
  static size_t intervalsPerYear;
  // 1.0/intervalsPerYear
  static double yearsPerInterval;
  // Maximum age of individuals in a scenario in time intervals
  static int maxAgeIntervals;
  //@}
  
  ///@brief Command-line options
  //@{
  static CLO::CLO clOptions;
  static string clResourcePath;
  static bool compressCheckpoints;
  //@}
  
  /** @brief Checkpointing.
   *
   * Not really required; more to confirm things are expected. */
  void read (istream& in);
  void write (ostream& out);	///< ditto
  
  /** Prepend if path is relative, prepend it with clResourcePath.
   * Then passes the resulting (or original) path through
   * BoincWrapper::resolveFile() and returns the result. */
  static string lookupResource (const string& path);
  
private:
  /// Sets modelVersion, checking for incompatible versions.
  static void setModelVersion ();
};

#endif
