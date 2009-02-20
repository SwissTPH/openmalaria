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

#ifndef Hsimulation
#define Hsimulation
#include <fstream>

extern int simulationDuration;

/*
  The estimated proportion of total simulation time we spend in mainSimulation.
  Used to give BOINC a hint on the progress.
*/
extern double relTimeInMainSim;

class TransmissionModel;
class Population;
class Summary;

static const char* const checkpoint_name = "checkpoint";
static const char* const checkpoint_0_name = "checkpoint0";
static const char* const checkpoint_1_name = "checkpoint1";

//! Main simulation class
class Simulation{

 public: 
 
  //Time since start of the simulation (TODO: Temporary step in refactoring. This variable should not be public"
  static int simulationTime;
  
  //Time counter during the main simulation. 
  static int timeStep;

  // Summary generator 
  static Summary* gMainSummary;

  //!  Inititalise all step specific constants and variables.
  Simulation();
  ~Simulation();

  //! Entry point to simulation.
  int start();

  //!  This procedure starts with the current state of the simulation 
  /*! It continues updating    assuming:
    (i)		the default (exponential) demographic model
    (ii)	the entomological input defined by the EIRs in intEIR()
    (iii)	the intervention packages defined in Intervention()
    (iv)	the survey times defined in Survey() */
  void mainSimulation();

  /*! Run the simulation using the equilibrium inoculation rates over one complete
    lifespan (maxAgeIntervals) to reach immunological equilibrium in all age
    classes. Don't report any events */
  void updateOneLifespan();

 private:

  TransmissionModel* _transmissionModel;
  Population* _population;

  //! This is the checkpoint data (the snapshot)
  std::fstream f_checkpoint;

  void validateInput();
  void tearDown();

   //!  Write the current checkpoint to a file.
  void writeCheckpoint();

  //! Open the checkpoint file
  void readCheckpoint();

  //! This function reads the checkpoint from the file in which is written
  /*! if it exists and initializes all the data return true if the workunit is
    checkpointed and false otherwise */
  bool isCheckpoint();

};

#endif
