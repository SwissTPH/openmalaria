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
#include "util/random.h"

#include <fstream>
#include "gzstream.h"


namespace OM {

// -----  Set-up & tear-down  -----

Simulation::Simulation(util::Checksum ck)
: simPeriodEnd(0), totalSimDuration(0), phase(STARTING_PHASE), _population(NULL), workUnitIdentifier(0), cksum(ck)
{
    Global::simulationTime = 0;
    Global::timeStep = numeric_limits<int>::min();
    
    // Initialize input variables and allocate memory.
    // We try to make initialization hierarchical (i.e. most classes initialise
    // through Population::init).
    gsl::setUp (InputData.getISeed());
    util::random::seed (InputData.getISeed());
    util::ModelOptions::init ();
    Surveys.init();
    Population::init();
    _population = new Population();
    
    workUnitIdentifier = InputData.getScenario().getWuID();
}

Simulation::~Simulation(){
  //free memory
  Population::clear();
  
  gsl::tearDown();
}


// -----  run simulations  -----

int Simulation::start(){
    // +1 to let final survey run
    totalSimDuration = _population->_transmissionModel->vectorInitDuration() + Global::maxAgeIntervals + Surveys.getFinalTimestep() + 1;
    int testCheckpointStep = -1;	// declare here so goto doesn't cross initialization
    
    if (isCheckpoint()) {
	readCheckpoint();
    } else {
	_population->createInitialHumans();
    }
    
    
    // phase loop
    while (true) {
	// loop for steps within a phase
	while (Global::simulationTime < simPeriodEnd) {
	    // checkpoint
	    if (util::BoincWrapper::timeToCheckpoint() || Global::simulationTime == testCheckpointStep) {
		writeCheckpoint();
		util::BoincWrapper::checkpointCompleted();
		if (Global::simulationTime == testCheckpointStep)
		    throw util::cmd_exit ("Checkpoint test: checkpoint written");
	    }
	    
	    // interventions
	    if (Global::timeStep == Surveys.currentTimestep) {
		_population->newSurvey();
		Surveys.incrementSurveyPeriod();
	    }
	    _population->implementIntervention(Global::timeStep);
	    
	    // update
	    ++Global::simulationTime;
	    _population->update1();
	    ++Global::timeStep;
	    util::BoincWrapper::reportProgress (double(Global::simulationTime) / totalSimDuration);
	}
	
	// phase transition: end of one phase to start of next
	if (phase == STARTING_PHASE) {
	    simPeriodEnd = _population->_transmissionModel->vectorInitDuration();
	    phase = VECTOR_FITTING;
	} else if (phase == VECTOR_FITTING) {
	    int extend = _population->_transmissionModel->vectorInitIterate ();
	    if (extend > 0) {	// repeat phase
		simPeriodEnd += extend;
		totalSimDuration += extend;
	    } else {		// next phase
		simPeriodEnd += Global::maxAgeIntervals;
		phase = ONE_LIFE_SPAN;
	    }
	} else if (phase == ONE_LIFE_SPAN) {
	    simPeriodEnd = totalSimDuration;
	    Global::timeStep=0;
	    _population->preMainSimInit();
	    _population->newSurvey();	// Only to reset TransmissionModel::innoculationsPerAgeGroup
	    Surveys.incrementSurveyPeriod();
	    phase = MAIN_PHASE;
	} else if (phase == MAIN_PHASE) {
	    phase = END_SIM;
	    break;
	}
	testCheckpointStep = -1;
	if (util::CommandLine::option (util::CommandLine::TEST_CHECKPOINTING))
	    testCheckpointStep = Global::simulationTime + (simPeriodEnd - Global::simulationTime) / 2;
    }
    
    delete _population;	// must destroy all Human instances to make sure they reported past events
    Surveys.writeSummaryArrays();
    
    // Write scenario checksum, only if simulation completed.
    cksum.writeToFile (util::BoincWrapper::resolveFile ("scenario.sum"));
    
    return 0;
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
    util::checkpoint::header (stream);
    util::CommandLine::staticCheckpoint (stream);
    util::random::checkpoint (stream);
    Population::staticCheckpoint (stream);
    Surveys & stream;
    
    Global::simulationTime & stream;
    Global::timeStep & stream;
    simPeriodEnd & stream;
    totalSimDuration & stream;
    phase & stream;
    (*_population) & stream;
    
    // Check scenario.xml and checkpoint files correspond:
    int oldWUID(workUnitIdentifier);
    util::Checksum oldCksum(cksum);
    workUnitIdentifier & stream;
    cksum & stream;
    if (workUnitIdentifier != oldWUID || cksum != oldCksum)
	throw util::checkpoint_error ("mismatched checkpoint");
    
    stream.ignore (numeric_limits<streamsize>::max()-1);	// skip to end of file
    if (stream.gcount () != 0) {
	ostringstream msg;
	msg << "Checkpointing file has " << stream.gcount() << " bytes remaining." << endl;
	throw util::checkpoint_error (msg.str());
    } else if (stream.fail())
	throw util::checkpoint_error ("stream read error");
}

void Simulation::checkpoint (ostream& stream) {
    util::checkpoint::header (stream);
    if (stream == NULL || !stream.good())
	throw util::checkpoint_error ("Unable to write to file");
    util::timer::startCheckpoint ();
    
    util::CommandLine::staticCheckpoint (stream);
    util::random::checkpoint (stream);
    Population::staticCheckpoint (stream);
    Surveys & stream;
    
    Global::simulationTime & stream;
    Global::timeStep & stream;
    simPeriodEnd & stream;
    totalSimDuration & stream;
    phase & stream;
    (*_population) & stream;
    workUnitIdentifier & stream;
    cksum & stream;
    
    util::timer::stopCheckpoint ();
    if (stream.fail())
	throw util::checkpoint_error ("stream write error");
}

}