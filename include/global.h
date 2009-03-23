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
#ifndef Hmod_global
#define Hmod_global
#include "Constant.h"
#include <fcntl.h>
#include <math.h>
#include <vector>
using namespace std;

class Global
{
public:
  /// Sets parameters in Global.
  static void initGlobal ();
  
  static int modIntervalsPerYear (int i);

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
   //Simulation time steps per year
  static int intervalsPerYear;
   //Maximum age of individuals in a scenario in time intervals
  static int maxAgeIntervals;
  static int simulationMode;
   //pre-erythrocytic latent period, in time steps
  static int latentp;
  
/*
  Size of the human population
  Moved from population.f so that transmission model.f can see it.
*/
  static vector<int> infantDeaths;
  static vector<int> infantIntervalsAtRisk;
  //@}
  
private:
  /// Sets modelVersion, checking for incompatible versions.
  static void setModelVersion ();
};

inline int isOptionIncluded (int allOptions, int option) {
  return allOptions & (1 << option);
};

#endif
