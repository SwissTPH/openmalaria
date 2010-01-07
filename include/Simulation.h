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

#ifndef Hmod_Simulation
#define Hmod_Simulation

#include "Global.h"
#include "Population.h"

using namespace std;

namespace OM {
    
//! Main simulation class
class Simulation{
public: 
    //!  Inititalise all step specific constants and variables.
    Simulation(util::Checksum);
    ~Simulation();
    
    //! Entry point to simulation.
    int start();
    
private:
    /** @brief checkpointing functions
    *
    * readCheckpoint/writeCheckpoint prepare to read/write the file,
    * and read/write read and write the actual data. */
    //@{
    //! This function reads the checkpoint from the file in which is written
    /*! if it exists and initializes all the data return true if the workunit is
    checkpointed and false otherwise */
    bool isCheckpoint();
    
    void writeCheckpoint();
    void readCheckpoint();
    
    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);
    //@}
    
    
    // Data
    int simPeriodEnd;
    int totalSimDuration;
    enum Phase {
	STARTING_PHASE = 0,
	/// Initialisation phase for Vector model
	VECTOR_FITTING,
	/*! Run the simulation using the equilibrium inoculation rates over one complete
	    lifespan (maxAgeIntervals) to reach immunological equilibrium in all age
	    classes. Don't report any events */
	ONE_LIFE_SPAN,
	//!  This procedure starts with the current state of the simulation 
	/*! It continues updating    assuming:
	    (i)		the default (exponential) demographic model
	    (ii)	the entomological input defined by the EIRs in intEIR()
	    (iii)	the intervention packages defined in Intervention()
	    (iv)	the survey times defined in Survey() */
	MAIN_PHASE,
	END_SIM		// should have largest value of all enumerations
    };
    int phase;
    
    Population* _population;
    
    /** Some identifier is needed to prevent checkpoint cheats. Ideally a unique identifier per
     * workunit, but a random integer number should do the job. */
    int workUnitIdentifier;
    
    /// Checksum of scenario; also used for work-unit validation on checkpoint loading.
    util::Checksum cksum;
    
    friend class VectorAnophelesSuite;
};

}
#endif
