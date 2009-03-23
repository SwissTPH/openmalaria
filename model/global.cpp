/*

 This file is part of OpenMalaria.
 
 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include <cstdlib>
#include "global.h"
#include "inputData.h"
#include <cmath>
#include <iostream>

/*
Contains global variables and constants and utility functions that are used in different modules.

Constants (parameter)

*/

ModelVersion Global::modelVersion;
int Global::interval;
int Global::intervalsPerYear;
int Global::maxAgeIntervals;
int Global::simulationMode;
int Global::latentp;

vector<int> Global::infantIntervalsAtRisk;
vector<int> Global::infantDeaths;

void Global::setModelVersion () {
  modelVersion = (ModelVersion) get_model_version();
  /* To print flags as binary:
  cerr << "Model version: ";
  for (int i = 24; i >= 0; --i)
    cerr << ((modelVersion >> i) & 1);
  cerr << endl; */
  
  // Or'd flags of incompatibility triangle from
  // "description of variables for interface" excel sheet
  const int INCOMPATIBLITITIES[NUM_VERSIONS] = {
    0,
    WITHIN_HOST_PARASITE | INCLUDES_PK_PD,	// 1
    LOGNORMAL_MASS_ACTION | LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM,	// 2
    WITHIN_HOST_PARASITE | INCLUDES_PK_PD | TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET,	// 3
    LOGNORMAL_MASS_ACTION_PLUS_PRE_IMM,	// 4
    TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET,	// 5
    WITHIN_HOST_PARASITE | INCLUDES_PK_PD | TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET,	// 6
    WITHIN_HOST_PARASITE | INCLUDES_PK_PD,	// 7
    WITHIN_HOST_PARASITE,	// 8
    0,	// 9
    0,	// 10
    MUELLER_MORBIDITY_MODEL,	// 11
    0,	// 12
    0,	// 13
    0,	// 14
    COMORB_TRANS_HET | TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,	// 15
    COMORB_TRANS_HET | TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,	// 16
    COMORB_TRANS_HET | TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,	// 17
    TRANS_TREAT_HET | COMORB_TREAT_HET | TRIPLE_HET,	// 18
    COMORB_TREAT_HET | TRIPLE_HET,	// 19
    TRIPLE_HET,	// 20
    0,	// 21
    0,
    0
  };
  
  for (size_t i = 0; i < NUM_VERSIONS; ++i)
    if (((modelVersion >> i) & 1) &&
          modelVersion & INCOMPATIBLITITIES[i]) {
      cerr << "Incompatible model versions: flag " << i << " is incompatible with other flags." << endl;
      throw 0;
    }
  if (modelVersion & (ATTENUATION_ASEXUAL_DENSITY | MAX_DENS_CORRECTION | INNATE_MAX_DENS | MAX_DENS_RESET))
    cerr << "Warning: model version used is deprecated" << endl;
}

void Global::initGlobal () {
  setModelVersion();
  interval=get_interval();
  if (daysInYear % interval !=  0) {
    cerr << "daysInYear not a multiple of interval" << endl;
    exit(-12);
  }
  intervalsPerYear = daysInYear/interval;
  infantDeaths.resize(intervalsPerYear);
  infantIntervalsAtRisk.resize(intervalsPerYear);
  latentp=get_latentp();
  maxAgeIntervals=(int)get_maximum_ageyrs()*intervalsPerYear;
}

int Global::modIntervalsPerYear (int i) {
    int valmodIntervalsPerYear = i % intervalsPerYear;
    if ( valmodIntervalsPerYear ==  0) {
        valmodIntervalsPerYear=intervalsPerYear;
    }
    return valmodIntervalsPerYear;
}
