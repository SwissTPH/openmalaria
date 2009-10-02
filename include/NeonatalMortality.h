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

#ifndef Hmod_NeonatalMortality
#define Hmod_NeonatalMortality

#include "global.h"
#include <list>

class Human;

class NeonatalMortality {
public:
  /// Initialisation
  static void init();
  
  /** Called for each birth; returns true if infant dies due to mother's
   * infection. */
  static bool eventNeonatalMortality();
  
  /** Calculate risk of a neonatal mortality based on humans 20-25 years old. */
  static void update (const list<Human>& population);
  
private:
  /** Calculates the risk of neonatal mortality. */
  static void calculateRiskFromMaternalInfection(int nCounter, int pCounter);
  
  /** Probability for a newborn to die (indirect death) because the mother is
   * infected. Depends on the prevalence of parasitaemia in mother at some
   * previous t. */
  static double _riskFromMaternalInfection;
  //! array for stored prevalences 20-25 years for 5 months (for neonatal deaths)
  static std::vector<double> _prevalenceByGestationalAge;
};

#endif
