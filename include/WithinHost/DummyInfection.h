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

#ifndef Hmod_DummyInfection
#define Hmod_DummyInfection
#include "WithinHost/Infection.h"
#include <fstream>
#include <fcntl.h>

//!  Models of infection.
/*!
  Models related to the within-host dynamics of infections.
*/
class DummyInfection : public Infection {
public:
  //! Constructor
  DummyInfection(int simulationTime);
  
  /** Checkpoint-reading constructor */
  DummyInfection (istream& in);
  
  /** Destructor
   * 
   * NOTE: this destructor does nothing to allow shallow copying to the
   * population list. destroy() does the real freeing and must be
   * called explicitly. */
  ~DummyInfection();
  void destroy();
  
  void write (ostream& out) const;
  
  static void init ();
  
  //! Get the last timestep before the infection is cleared.
  /*!
    \return The interval before clearance.
  */
  int getEndDate();

  /// Multiplies the density by x.
  void multiplyDensity(double x) { _density *= x; };
  //! Get the density of the infection
  double getDensity() { return _density; };

  //! Start date of the infection
  int getStartDate() { return _startdate; };
  
  /*
  //! Dummy function. Should override \sa determineDensities.
  */ 
  void determineWithinHostDensity();
  
private:
  //! Arbitrary maximum duration of the infection
  int _duration; 
};

#endif

