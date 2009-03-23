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

#ifndef Hmod_morbidity
#define Hmod_morbidity

#include <iostream>

using namespace std;

/*! Morbidity Model abstract base class. */
class MorbidityModel {
public:
  virtual double getPEpisode(double timeStepMaxDensity, double totalDensity)=0;
  virtual double getPyrogenThres();
  virtual void write(ostream& out) const=0;
  virtual void read(istream& in)=0;
  
  // Static:
  /// Calls static init on all MorbidityModels.
  static void initModels();
  static MorbidityModel* createMorbidityModel();
};

#endif
