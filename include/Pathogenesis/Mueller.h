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

#ifndef Hmod_muellerMorb
#define Hmod_muellerMorb

#include <iostream>
#include "global.h"
#include "Pathogenesis/PathogenesisModel.h"

using namespace std;

/*! Mueller presentation model.
*/
class MuellerPathogenesis : public PathogenesisModel {
public:
  MuellerPathogenesis(double cF) :
    PathogenesisModel(cF) {}
  MuellerPathogenesis(istream& in) :
    PathogenesisModel(in) {}
  
  virtual double getPEpisode(double timeStepMaxDensity, double totalDensity);

  // Static:
  static void init();

private:
  static double rateMultiplier_31;
  static double densityExponent_32;
};

#endif
