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

#include "Simulation.h"

#include "util/BoincWrapper.h"
#include "util/timer.h"
#include "util/gsl.h"
#include "Surveys.h"
#include "Global.h"
#include "Transmission/TransmissionModel.h"
#include "inputData.h"
#include "util/CommandLine.hpp"
#include "util/ModelOptions.hpp"
#include "util/errors.hpp"

#include <fstream>
#include "gzstream.h"


namespace OM {

// -----  Set-up & tear-down  -----

Simulation::Simulation()
{
    Global::simulationTime = 0;
    Global::timeStep = -1;
    
    // Initialize input variables and allocate memory.
  // We try to make initialization hierarchical (i.e. most classes initialise
  // through Population::init).
  gsl::setUp();
  util::ModelOptions::set (InputData.get_model_version());
  Surveys.init();
  Population::init();
  _population = new Population();
}

Simulation::~Simulation(){
  //free memory
  Population::clear();
  
  gsl::tearDown();
}


// -----  run simulations  -----

int Simulation::start(){
  _population->estimateRemovalRates();
  if (isCheckpoint()) {
    _population->setupPyramid(true);
    readCheckpoint();
  } else {
    _population->setupPyramid(false);
  }
  
  simPeriodEnd = _population->_transmissionModel->vectorInitDuration();
  // +1 to let final survey run
  totalSimDuration = simPeriodEnd + Global::maxAgeIntervals + Surveys.getFinalTimestep() + 1;
  
  while (Global::simulationTime < simPeriodEnd) {
    vectorInitialisation();
    int extend = _population->_transmissionModel->vectorInitIterate ();
    simPeriodEnd += extend;
    totalSimDuration += extend;
  }
  
  simPeriodEnd += Global::maxAgeIntervals;
  updateOneLifespan();
  
  simPeriodEnd = totalSimDuration;
  mainSimulation();
  return 0;
}

void Simulation::vectorInitialisation () {
  while(Global::simulationTime < simPeriodEnd) {
    ++Global::simulationTime;
    _population->update1();
    
    util::BoincWrapper::reportProgress (double(Global::simulationTime) / totalSimDuration);
    if (util::BoincWrapper::timeToCheckpoint()) {
      writeCheckpoint();
      util::BoincWrapper::checkpointCompleted();
    }
  }
}

void Simulation::updateOneLifespan () {
  int testCheckpointStep = -1;
  if (util::CommandLine::option (util::CommandLine::TEST_CHECKPOINTING))
    testCheckpointStep = simPeriodEnd - Global::maxAgeIntervals / 2;
  while(Global::simulationTime < simPeriodEnd) {
    ++Global::simulationTime;
    _population->update1();
    
    util::BoincWrapper::reportProgress (double(Global::simulationTime) / totalSimDuration);
    if (util::BoincWrapper::timeToCheckpoint() || Global::simulationTime == testCheckpointStep) {
      writeCheckpoint();
      util::BoincWrapper::checkpointCompleted();
      if (util::CommandLine::option (util::CommandLine::TEST_CHECKPOINTING))
	throw util::cmd_exit ("Checkpoint test: written checkpoint");
    }
  }
}

void Simulation::mainSimulation(){
  Global::timeStep=0;
  _population->preMainSimInit();
  _population->newSurvey();	// Only to reset TransmissionModel::innoculationsPerAgeGroup
  Surveys.incrementSurveyPeriod();
  
  while(Global::simulationTime < simPeriodEnd) {
    if (Global::timeStep == Surveys.currentTimestep) {
      _population->newSurvey();
      Surveys.incrementSurveyPeriod();
    }
    _population->implementIntervention(Global::timeStep);
    //Calculate the current progress
    util::BoincWrapper::reportProgress(double(Global::simulationTime) / totalSimDuration);
    
    ++Global::simulationTime;
    _population->update1();
    ++Global::timeStep;
    //FIXME: also write checkpoints here.
  }
  delete _population;	// must destroy all Human instances to make sure they reported past events
  Surveys.writeSummaryArrays();
}


// -----  checkpointing: set up read/write stream  -----

const char* CHECKPOINT = "checkpoint";

bool Simulation::isCheckpoint(){
  ifstream checkpointFile(CHECKPOINT,ios::in);
  // If not open, file doesn't exist (or is inaccessible)
  return checkpointFile.is_open();
}
int readCheckpointNum () {
    ifstream checkpointFile;
    checkpointFile.open(CHECKPOINT, fstream::in);
    int checkpointNum;
    checkpointFile >> checkpointNum;
    checkpointFile.close();
    if (!checkpointFile)
	throw util::checkpoint_error ("error reading from file \"checkpoint\"");
    return checkpointNum;
}

void Simulation::writeCheckpoint(){
  // We alternate between two checkpoints, in case program is closed while writing.
  const int NUM_CHECKPOINTS = 2;
  
  int checkpointNum = 0;
    if (isCheckpoint()) {
	checkpointNum = readCheckpointNum();
	// Get next checkpoint number:
	checkpointNum = (checkpointNum + 1) % NUM_CHECKPOINTS;
    }
  
  gsl::rngSaveState (checkpointNum);
  
  // Open the next checkpoint file for writing:
  ostringstream name;
  name << CHECKPOINT << checkpointNum;
  if (util::CommandLine::option (util::CommandLine::COMPRESS_CHECKPOINTS)) {
    name << ".gz";
    ogzstream out(name.str().c_str(), ios::out | ios::binary);
    checkpoint (out);
    out.close();
  } else {
    ofstream out(name.str().c_str(), ios::out | ios::binary);
    checkpoint (out);
    out.close();
  }
  
  {	// Indicate which is the latest checkpoint file.
    ofstream checkpointFile;
    checkpointFile.open(CHECKPOINT,ios::out);
    checkpointFile << checkpointNum;
    checkpointFile.close();
    if (!checkpointFile)
	throw util::checkpoint_error ("error writing to file \"checkpoint\"");
  }
}

void Simulation::readCheckpoint() {
    int checkpointNum = readCheckpointNum();
    
  // Open the latest file
  ostringstream name;
  name << CHECKPOINT << checkpointNum;	// try uncompressed
  ifstream in(name.str().c_str(), ios::in | ios::binary);
  if (in.good()) {
    checkpoint (in);
    in.close();
  } else {
    name << ".gz";				// then compressed
    igzstream in(name.str().c_str(), ios::in | ios::binary);
    if (!in.good())
      throw util::checkpoint_error ("Unable to read file");
    checkpoint (in);
    in.close();
  }
  
  gsl::rngLoadState (checkpointNum);
  cerr << "Loaded checkpoint from: " << name.str() << endl;
  
  // On resume, write a checkpoint so we can tell whether we have identical checkpointed state
  if (util::CommandLine::option (util::CommandLine::TEST_CHECKPOINTING))
    writeCheckpoint();
}


//   -----  checkpointing: Simulation data  -----

void Simulation::checkpoint (istream& stream) {
    util::CommandLine::staticCheckpoint (stream);
    // FIXME: appears problematic (though may not be needed):
    //   Population::staticRead(stream);
    Surveys & stream;
    
    Global::simulationTime & stream;
    Global::timeStep & stream;
    simPeriodEnd & stream;
    totalSimDuration & stream;
    (*_population) & stream;
    
    // Read trailing white-space (final endl has not yet been read):
    while (stream.good() && isspace (stream.peek()))
	stream.get();
    if (!stream.eof()) {	// if anything else is left
	cerr << "Error (checkpointing): not the whole checkpointing file was read;";
	ifstream *ifCP = dynamic_cast<ifstream*> (&stream);
	if (ifCP) {
	    streampos i = ifCP->tellg();
	    ifCP->seekg(0, ios_base::end);
	    cerr << ifCP->tellg()-i << " bytes remaining:";
	    ifCP->seekg (i);
	} else	// igzstream can't seek
	    cerr << " remainder:" << endl;
	cerr << endl << stream.rdbuf() << endl;
    }
}

void Simulation::checkpoint (ostream& stream) {
    if (stream == NULL || !stream.good())
	throw util::checkpoint_error ("Unable to write to file");
    timer::startCheckpoint ();
    stream.precision(20);
    
    util::CommandLine::staticCheckpoint (stream);
    //   Population::staticWrite(stream);
    Surveys & stream;
    
    Global::simulationTime & stream;
    Global::timeStep & stream;
    simPeriodEnd & stream;
    totalSimDuration & stream;
    (*_population) & stream;
    
    timer::stopCheckpoint ();
    if (stream.fail())
	throw util::checkpoint_error ("stream write error");
}

}