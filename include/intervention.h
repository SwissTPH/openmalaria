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
#ifndef Hmod_intervention
#define Hmod_intervention
#include <fcntl.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"

// Forward declare: doesn't need to be known about here
namespace scnXml {
  class VaccineDescription;
}

/** Vaccine intervention parameters.
 * 
 * Used to represent PEV, BSV and TBV vaccines. */
class Vaccine {
public:
  // Static:
  /// Set parameters from xml
  static void initParameters ();
  /// Free memory
  static void clearParameters ();
  
  /// True if any types of vaccine are in use.
  static bool anyVaccine;
  
  /*! Common to all vaccine types. Number of vaccine doses that are given
   * either through EPI or as EPI Boosters. */
  static size_t _numberOfEpiDoses;
  
   /** Target age for EPI-like vaccination, in time steps. */
  static vector<int> targetAgeTStep;
   ///Coverage , as a proportion of the poulation in the target age range
  static double *vaccineCoverage;
  
  /// Preerythrocytic reduces h vaccine parameters
  static Vaccine PEV;
  /// Erythrocytic reduces y vaccine parameters
  static Vaccine BSV;
  /// Transmission blocking reduces k vaccine parameters
  static Vaccine TBV;
  
  //Non-static:
  Vaccine() : active(false), decay(1.0) {}
  
  /// True if this vaccine is in use
  bool active;
  
  /** Get the efficacy of the vaccine.
   * 
   * @param numPrevDoses The number of prior vaccinations of the individual. */
  double getEfficacy (int numPrevDoses);
  
  /// exp(-Decay rate)
  double decay;
private:
  /** Per-type initialization
   * @returns decay */
  void initVaccine (const scnXml::VaccineDescription* vd);
  
  /* Vaccine type specific parameters
   * Initial mean efficacy, definition depends on vaccine type */
  vector<double> initialMeanEfficacy;
  // Distribution of efficacies among individuals, parameter to sample from beta dist.
  double efficacyB;
};

#endif
