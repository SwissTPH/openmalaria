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
#include <fcntl.h>
#include "intervention.h"
#include "WithinHost/Infection.h"
#include <fstream>

//!  Models of infection.
/*!
  Models related to the within-host dynamics of infections.
*/
class DummyInfection : public Infection {
 //TODO: should be private, and immune decay and immune proxies need to be discussed in light of new WIH-models 
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
  
  int getDuration() { return _duration; };

  /*
  //! Dummy function. Should override \sa determineDensities.
  */ 
  void determineWithinHostDensity();

  //! Initialises infection duration.
  /*! 
    Initialises infection duration sampling from log normal distribution using parameters for 53 patients from Georgia.
    Mean log duration of an infection values from AJTM p.9 eq.5.
    \return The duration in simulation intervals.
  */
  int infectionDuration();

  
  //! Write an infection to the checkpointing file.
  /*!
    \param funit Checkpoint file.
  */
  void writeInfectionToFile (fstream& funit);

  float getCumulativeHstar() const {return cumulativeHstar;};
  float getCumulativeYstar() const {return cumulativeYstar;};
};

#endif

