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

#ifndef Hmod_whost
#define Hmod_whost

#include "global.h"
#include <iostream>

using namespace std;

class Human;

/*! Within Host Model abstract class.
 * Dont forget to create friend << and >> for subclasses.
 */
class WithinHostModel {
  protected:
    virtual void write(ostream& out) const =0;
    virtual void read(istream& in) =0;

  public:
    WithinHostModel() {}

    friend ostream& operator<<(ostream& out, const WithinHostModel &model);
    friend istream& operator>>(istream& in, WithinHostModel &model);


    virtual void calculateDensities(Human&) =0;

};

#endif
