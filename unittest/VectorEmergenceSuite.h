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

#ifndef Hmod_VectorEmergenceSuite
#define Hmod_VectorEmergenceSuite

#include <cxxtest/TestSuite.h>
#include "Transmission/VectorEmergence.h"
#include "Transmission/VectorSpecies.h"

// We want to hide normal output, so route it here instead of cout
ofstream null("\0");

class VectorEmergenceSuite : public CxxTest::TestSuite
{
public:
  VectorEmergenceSuite () :
      emerge(3, 10, POP_SIZE, AVG_AVAIL, 1.6, 0.33, .95, .95, .94, .93,
	     YEAR_LEN, null, "\0")
  {
  }
  
  void testDummy () {
    TS_WARN ("TODO: write VectorEmergence tests");
  }
  
  void testWholeCalculation () {
    Global::clOptions = static_cast<CLO::CLO> (Global::clOptions | CLO::ENABLE_ERC);
    
    vector<double> humanInfectivityInit(YEAR_LEN, 0.0);	//TODO: init
    vector<double> EIRInit(YEAR_LEN, 0.0);
    // Fourier cooeficients for EIR array
    vector<double> fc (5,0.0);
    fc[0] = -0.926517;	// a0
    fc[1] = -0.692164;	fc[2] =  0.002098;	// a1,b1
    fc[3] =  0.401189;	fc[4] = -0.375356;	// a2,b2
    VectorTransmissionSpecies::calcInverseDFTExp (EIRInit, fc);
    
    double emergenceRate[YEAR_LEN];
    const double temp = POP_SIZE*POP_SIZE*AVG_AVAIL;
    for (int i = 0; i < YEAR_LEN; i++) {
      emergenceRate[i] = EIRInit[i]*temp;
    }
    
    emerge.CalcInitMosqEmergeRate(1,1,	//NOTE: no support for these not-being 1 yet
				  &humanInfectivityInit[0],
				  EIRInit,
				  emergenceRate);
    TS_WARN ("TODO: test output");
  }
  
private:
  // Number of "days" in our "year" (to speed up tests)
  static const int YEAR_LEN = 10;
  // Population size
  static const int POP_SIZE = 1000;
  // Average availability
  static const double AVG_AVAIL = 0.0072;
  
  VectorEmergence emerge;
};

#endif
