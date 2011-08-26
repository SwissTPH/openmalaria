/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "Simulation.h"

#include "util/BoincWrapper.h"
#include "util/timer.h"
#include "Monitoring/Continuous.h"
#include "Monitoring/Surveys.h"
#include "Global.h"
#include "Transmission/TransmissionModel.h"
#include "inputData.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "PopulationStats.h"

#include <fstream>
#include <gzstream/gzstream.h>


namespace OM {
    using Monitoring::Continuous;
    using Monitoring::Surveys;

// -----  Set-up & tear-down  -----

Simulation::Simulation(util::Checksum ck) :
    simPeriodEnd(0),
    totalSimDuration(0),
    phase(STARTING_PHASE),
    workUnitIdentifier(0),
    cksum(ck)
{
    OM::TimeStep::init(InputData().getModel().getParameters().getInterval(),
                       InputData().getDemography().getMaximumAgeYrs());
    
    // Initialize input variables and allocate memory.
    // We try to make initialization hierarchical (i.e. most classes initialise
    // through Population::init).
    util::random::seed (InputData().getModel().getParameters().getIseed());
    util::ModelOptions::init ();
    Surveys.init();
    Population::init();
    population = auto_ptr<Population>(new Population);
    interventions = auto_ptr<InterventionManager>( new InterventionManager( InputData().getInterventions(), *population ) );
    Surveys.initCohortOnly( *interventions );
    
    workUnitIdentifier = InputData().getWuID();
}

Simulation::~Simulation(){
    //free memory
    Population::clear();
}


// -----  run simulations  -----
enum Phase {
    STARTING_PHASE = 0,
    /*! Run the simulation using the equilibrium inoculation rates over one complete
        lifespan (maxAgeIntervals) to reach immunological equilibrium in all age
        classes. Don't report any events */
    ONE_LIFE_SPAN,
    /** Initialisation/fitting phase for transmission models. */
    TRANSMISSION_INIT,
    //!  This procedure starts with the current state of the simulation 
    /*! It continues updating    assuming:
        (i)         the default (exponential) demographic model
        (ii)        the entomological input defined by the EIRs in intEIR()
        (iii)       the intervention packages defined in Intervention()
        (iv)        the survey times defined in Survey() */
    MAIN_PHASE,
    END_SIM         // should have largest value of all enumerations
};

int Simulation::start(){
    TimeStep::simulation = TimeStep( 0 );
    
    // Make sure warmup period is at least as long as a human lifespan, as the
    // length required by vector warmup, and is a whole number of years.
    TimeStep humanWarmupLength = TimeStep::maxAgeIntervals;
    if( humanWarmupLength < population->_transmissionModel->minPreinitDuration() ){
        cerr << "Warning: human life-span (" << humanWarmupLength.inYears();
        cerr << ") shorter than length of warm-up requested by" << endl;
        cerr << "transmission model (" << population->_transmissionModel->minPreinitDuration().inYears();
        cerr << "). Transmission may be unstable; perhaps use forced" << endl;
        cerr << "transmission (mode=2) or a longer life-span." << endl;
        humanWarmupLength = population->_transmissionModel->minPreinitDuration();
    }
    humanWarmupLength = TimeStep::fromYears((humanWarmupLength.asInt()-1) / TimeStep::stepsPerYear + 1);
    
    totalSimDuration = humanWarmupLength  // ONE_LIFE_SPAN
        + population->_transmissionModel->expectedInitDuration()
        + Surveys.getFinalTimestep() + TimeStep( 1 );   // MAIN_PHASE: surveys; +1 to let final survey run
    
    if (isCheckpoint()) {
        Continuous::init( true );
        readCheckpoint();
    } else {
        Continuous::init( false );
        population->createInitialHumans();
    }
    // Set to either a checkpointing timestep or min int value. We only need to
    // set once, since we exit after a checkpoint triggered this way.
    TimeStep testCheckpointStep = util::CommandLine::getNextCheckpointTime( TimeStep::simulation );
    TimeStep testCheckpointDieStep = testCheckpointStep;        // kill program at same time
    
    // phase loop
    while (true) {
        // loop for steps within a phase
        while (TimeStep::simulation < simPeriodEnd) {
            // checkpoint
            if (util::BoincWrapper::timeToCheckpoint() || TimeStep::simulation == testCheckpointStep) {
                writeCheckpoint();
                util::BoincWrapper::checkpointCompleted();
            }
            if (TimeStep::simulation == testCheckpointDieStep)
                throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
            
            // do reporting (continuous and surveys)
            Continuous::report();
            if (TimeStep::interventionPeriod == Surveys.currentTimestep) {
                population->newSurvey();
                Surveys.incrementSurveyPeriod();
            }
            
            // deploy interventions
            interventions->deploy( *population );
            
            // update
            ++TimeStep::simulation;
            ++TimeStep::interventionPeriod;
            population->update1();
            
            util::BoincWrapper::reportProgress (double(TimeStep::simulation.asInt()) / totalSimDuration.asInt());
        }
        
        ++phase;        // advance to next phase
        if (phase == ONE_LIFE_SPAN) {
            simPeriodEnd = humanWarmupLength;
        } else if (phase == TRANSMISSION_INIT) {
            TimeStep iterate = population->_transmissionModel->initIterate();
            if( iterate > TimeStep(0) ) {
                simPeriodEnd += iterate;
                --phase;        // repeat phase
            } else {
                // nothing to do: start main phase immediately
            }
            // adjust estimation of final time step: end of current period + length of main phase
            totalSimDuration = simPeriodEnd + Surveys.getFinalTimestep() + TimeStep( 1 );
        } else if (phase == MAIN_PHASE) {
            // Start MAIN_PHASE:
            simPeriodEnd = totalSimDuration;
            TimeStep::interventionPeriod = TimeStep(0);
            population->preMainSimInit();
            population->newSurvey();       // Only to reset TransmissionModel::inoculationsPerAgeGroup
            Surveys.incrementSurveyPeriod();
        } else if (phase == END_SIM) {
            cerr << "sim end" << endl;
            break;
        }
        if (util::CommandLine::option (util::CommandLine::TEST_CHECKPOINTING)){
            // First of middle of next phase, or current value (from command line) triggers a checkpoint.
            TimeStep phase_mid = TimeStep::simulation + TimeStep( (simPeriodEnd - TimeStep::simulation).asInt() / 2 );
            // Don't checkpoint 0-length phases or do mid-phase checkpointing
            // when timed checkpoints were specified, and don't checkpoint
            // ONE_LIFE_SPAN phase if already past time humanWarmupLength:
            // these are extra transmission inits, and we don't want to
            // checkpoint every one of them.
            if( testCheckpointStep == TimeStep::never
                && phase_mid > TimeStep::simulation
                && (phase != ONE_LIFE_SPAN || TimeStep::simulation < humanWarmupLength)
            ){
                testCheckpointStep = phase_mid;
                // Test checkpoint: die a bit later than checkpoint for better
                // resume testing (specifically, ctsout.txt).
                testCheckpointDieStep = testCheckpointStep + TimeStep(2);
            }
        }
    }
    
    // Open a critical section; should prevent app kill while/after writing
    // output.txt, which we don't currently handle well.
    // Note: we don't end this critical section; we simply exit.
    util::BoincWrapper::beginCriticalSection();
    
    PopulationStats::print();
    
    population->flushReports();        // ensure all Human instances report past events
    Surveys.writeSummaryArrays();
    Continuous::finalise();
    
    // Write scenario checksum, only if simulation completed.
    cksum.writeToFile (util::BoincWrapper::resolveFile ("scenario.sum"));
    
# ifdef OM_STREAM_VALIDATOR
    util::StreamValidator.saveStream();
# endif
    
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
    int checkpointNum=0;
    checkpointFile >> checkpointNum;
    checkpointFile.close();
    if (!checkpointFile)
        throw util::checkpoint_error ("error reading from file \"checkpoint\"");
    return checkpointNum;
}

void Simulation::writeCheckpoint(){
    // We alternate between two checkpoints, in case program is closed while writing.
    const int NUM_CHECKPOINTS = 2;
    
    int oldCheckpointNum = 0, checkpointNum = 0;
    if (isCheckpoint()) {
        oldCheckpointNum = readCheckpointNum();
        // Get next checkpoint number:
        checkpointNum = (oldCheckpointNum + 1) % NUM_CHECKPOINTS;
    }
    
    {   // Open the next checkpoint file for writing:
        ostringstream name;
        name << CHECKPOINT << checkpointNum;
        //Writing checkpoint:
//      cerr << TimeStep::simulation << " WC: " << name.str();
        if (util::CommandLine::option (util::CommandLine::COMPRESS_CHECKPOINTS)) {
            name << ".gz";
            ogzstream out(name.str().c_str(), ios::out | ios::binary);
            checkpoint (out, checkpointNum);
            out.close();
        } else {
            ofstream out(name.str().c_str(), ios::out | ios::binary);
            checkpoint (out, checkpointNum);
            out.close();
        }
    }
    
    {   // Indicate which is the latest checkpoint file.
        ofstream checkpointFile;
        checkpointFile.open(CHECKPOINT,ios::out);
        checkpointFile << checkpointNum;
        checkpointFile.close();
        if (!checkpointFile)
            throw util::checkpoint_error ("error writing to file \"checkpoint\"");
    }
    // Truncate the old checkpoint to save disk space, when it existed
    if( oldCheckpointNum != checkpointNum
        && !util::CommandLine::option (
            util::CommandLine::TEST_DUPLICATE_CHECKPOINTS
        )       /* need original in this case */
    ) {
        ostringstream name;
        name << CHECKPOINT << oldCheckpointNum;
        if (util::CommandLine::option (util::CommandLine::COMPRESS_CHECKPOINTS)) {
            name << ".gz";
        }
        ofstream out(name.str().c_str(), ios::out | ios::binary);
        out.close();
    }
//     cerr << " OK" << endl;
}

void Simulation::readCheckpoint() {
    int checkpointNum = readCheckpointNum();
    
  // Open the latest file
  ostringstream name;
  name << CHECKPOINT << checkpointNum;  // try uncompressed
  ifstream in(name.str().c_str(), ios::in | ios::binary);
  if (in.good()) {
    checkpoint (in, checkpointNum);
    in.close();
  } else {
    name << ".gz";                              // then compressed
    igzstream in(name.str().c_str(), ios::in | ios::binary);
    //Note: gzstreams are considered "good" when file not open!
    if ( !( in.good() && in.rdbuf()->is_open() ) )
      throw util::checkpoint_error ("Unable to read file");
    checkpoint (in, checkpointNum);
    in.close();
  }
  
  // Keep size of stderr.txt minimal with a short message, since this is a common message:
  cerr <<TimeStep::simulation<<" RC"<<endl;
  
  // On resume, write a checkpoint so we can tell whether we have identical checkpointed state
  if (util::CommandLine::option (util::CommandLine::TEST_DUPLICATE_CHECKPOINTS))
    writeCheckpoint();
}


//   -----  checkpointing: Simulation data  -----

void Simulation::checkpoint (istream& stream, int checkpointNum) {
    try {
        util::checkpoint::header (stream);
        util::CommandLine::staticCheckpoint (stream);
        Population::staticCheckpoint (stream);
        Surveys & stream;
        Continuous::staticCheckpoint (stream);
#       ifdef OM_STREAM_VALIDATOR
        util::StreamValidator & stream;
#       endif
        
        TimeStep::interventionPeriod & stream;
        simPeriodEnd & stream;
        totalSimDuration & stream;
        phase & stream;
        (*population) & stream;
        PopulationStats::staticCheckpoint( stream );
        (*interventions) & stream;
        interventions->loadFromCheckpoint( *population, TimeStep::interventionPeriod );
        
        // read last, because other loads may use random numbers or expect time
        // to be negative
        TimeStep::simulation & stream;
        util::random::checkpoint (stream, checkpointNum);
        
        // Check scenario.xml and checkpoint files correspond:
        int oldWUID(workUnitIdentifier);
        util::Checksum oldCksum(cksum);
        workUnitIdentifier & stream;
        cksum & stream;
        if (workUnitIdentifier != oldWUID || cksum != oldCksum)
            throw util::checkpoint_error ("mismatched checkpoint");
    } catch (const util::checkpoint_error& e) { // append " (pos X of Y bytes)"
        ostringstream pos;
        pos<<" (pos "<<stream.tellg()<<" of ";
        stream.ignore (numeric_limits<streamsize>::max()-1);    // skip to end of file
        pos<<stream.tellg()<<" bytes)";
        throw util::checkpoint_error( e.what() + pos.str() );
    }
    
    
    stream.ignore (numeric_limits<streamsize>::max()-1);        // skip to end of file
    if (stream.gcount () != 0) {
        ostringstream msg;
        msg << "Checkpointing file has " << stream.gcount() << " bytes remaining." << endl;
        throw util::checkpoint_error (msg.str());
    } else if (stream.fail())
        throw util::checkpoint_error ("stream read error");
}

void Simulation::checkpoint (ostream& stream, int checkpointNum) {
    util::checkpoint::header (stream);
    if (stream == NULL || !stream.good())
        throw util::checkpoint_error ("Unable to write to file");
    util::timer::startCheckpoint ();
    
    util::CommandLine::staticCheckpoint (stream);
    Population::staticCheckpoint (stream);
    Surveys & stream;
    Continuous::staticCheckpoint (stream);
# ifdef OM_STREAM_VALIDATOR
    util::StreamValidator & stream;
# endif
    
    TimeStep::interventionPeriod & stream;
    simPeriodEnd & stream;
    totalSimDuration & stream;
    phase & stream;
    (*population) & stream;
    PopulationStats::staticCheckpoint( stream );
    (*interventions) & stream;
    
    TimeStep::simulation & stream;
    util::random::checkpoint (stream, checkpointNum);
    workUnitIdentifier & stream;
    cksum & stream;
    
    util::timer::stopCheckpoint ();
    if (stream.fail())
        throw util::checkpoint_error ("stream write error");
}

}
