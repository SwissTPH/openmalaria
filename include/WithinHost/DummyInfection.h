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

#include "WithinHost/CommonInfection.h"

namespace OM { namespace WithinHost {

//!  Models of infection.
/*!
  Models related to the within-host dynamics of infections.
*/
class DummyInfection : public CommonInfection {
public:
    /// For checkpointing (don't use for anything else)
    DummyInfection (istream& stream);
    //! Constructor
    DummyInfection (uint32_t protID);
    
    virtual ~DummyInfection () {}
  
  static void initParameters ();
  
  //! Get the last timestep before the infection is cleared.
  /*!
    \return The interval before clearance.
  */
  int getEndDate();

  
  virtual bool updateDensity (int simulationTime, double survivalFactor);
};

} }
#endif