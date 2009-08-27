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

#include "util/WeibullDecayedValue.h"
#include "inputData.h"

void WeibullDecayedValue::setParameters (double initial, double halflife, double k) {
  _initial = initial;
  _k = k;
  _constOverLambda = Global::yearsPerInterval / (halflife * pow(log(2.0), -1.0/k));
}
void WeibullDecayedValue::operator= (const scnXml::WeibullDecayedProportion& elt) {
  _initial = elt.getInitial ();
  _k = elt.getWeibullk().present() ?
    elt.getWeibullk().get() : 1.0;
  _constOverLambda = Global::yearsPerInterval / (elt.getHalflife() * pow(log(2.0), -1.0/_k));
}

double WeibullDecayedValue::operator() (int ageTSteps) {
  return _initial * exp(-pow(ageTSteps * _constOverLambda, _k));
}
