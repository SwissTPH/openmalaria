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

#ifndef Hmod_oldwhost
#define Hmod_oldwhost

#include <iostream>
#include "global.h"
#include "withinHostModel.h"

using namespace std;

class Human;
class Infection;

/*! Old Within Host Model class.
 */
class OldWithinHostModel : public WithinHostModel {
  private:
  //! time at which attenuated infection 'would' end if SP present
  int _SPattenuationt;
  double cumulativeY;
  double cumulativeh;
  double timeStepMaxDensity;
  void write(ostream& out) const;
  void read(istream& in);

  void SPAction();

  void calculateDensity(Infection *inf);

  public:
  OldWithinHostModel(Human *human);

  friend ostream& operator<<(ostream& out, const WithinHostModel &model);
  friend istream& operator>>(istream& in, WithinHostModel &model);


  void calculateDensities();

};

#endif
