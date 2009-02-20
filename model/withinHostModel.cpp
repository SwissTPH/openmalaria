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

#include "human.h"
#include "infection.h"
#include "withinHostModel.h"

using namespace std;


WithinHostModel::WithinHostModel(Human* human) {
  _human = human;
}

void WithinHostModel::calculateDensity(Infection *inf) {
}

void WithinHostModel::calculateDensities() {
}

void WithinHostModel::write(ostream& out) const {
}

void WithinHostModel::read(istream& in) {
}

ostream& operator<<(ostream& out, const WithinHostModel &model) {
  model.write(out);
  return out;
}

istream& operator>>(istream& in, WithinHostModel &model) {
  model.read(in);
  return in;
}

