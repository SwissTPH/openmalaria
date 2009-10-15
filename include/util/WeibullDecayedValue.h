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

#ifndef Hmod_WeibullDecayedValue
#define Hmod_WeibullDecayedValue

#include "Global.h"
#include "Simulation.h"

namespace scnXml {
  class WeibullDecayedProportion;
}

/** A class representing a value decayed by the Weibull cumulative density
 * function.
 * 
 * The age is not stored, thus if multiple objects need a property, each with
 * the same initial value and decay curves, they may share an instance of this
 * class, only passing the age to operator() to calculate their value.
 * 
 * The Weibull distribution can model exponential decay with k = 1.0, however
 * these calculations cannot be as fast as v2 = v1 * Const. */
class WeibullDecayedValue
{
public:
  /** Initialise values such that the value returned by operator() is always
   * zero when setParameters or operator= has not been called. */
  WeibullDecayedValue () :
    _initial(0.0), _k(1.0), _constOverLambda(0.0)
  {}
  
  /** Set the initial value and the Weibull distribution parameters.
   *
   * @param initial Initial value (returned by operator(0))
   * @param halflife Half-life for decay. λ = halflife * pow(ln(2), -1/k)
   * @param k Parameter k of the Weibull distribution; when 1.0 decay is
   *  exponential. */
  void setParameters (double initial, double halflife, double k = 1.0);
  /** As setParameters, but take values from the passed XML element. */
  void operator= (const scnXml::WeibullDecayedProportion&);
  
  /** Return the value decayed to age ageTSteps. */
  double operator() (int ageTSteps);
  
private:
  double _initial, _k;
  double _constOverLambda;	// Global::yearsPerInterval / lambda
};

#endif
