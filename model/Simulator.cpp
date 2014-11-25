/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#include "Simulator.h"

#include "Global.h"
#include "Transmission/TransmissionModel.h"
#include "PopulationStats.h"
#include "Parameters.h"
#include "Clinical/ClinicalModel.h"
#include "Monitoring/Continuous.h"
#include "interventions/InterventionManager.hpp"
#include "Population.h"
#include "WithinHost/Diagnostic.h"
#include "WithinHost/Genotypes.h"
#include "mon/management.h"
#include "util/BoincWrapper.h"
#include "util/timer.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "schema/scenario.h"

#include <fstream>
#include <gzstream/gzstream.h>


namespace OM {
    using Monitoring::Continuous;
    using interventions::InterventionManager;

bool Simulator::startedFromCheckpoint;  // static

const char* CHECKPOINT = "checkpoint";

enum Phase {
    STARTING_PHASE = 0,
    /** Run the simulation using the equilibrium inoculation rates over one
     * complete lifespan (sim::maxHumanAge()) to reach immunological
     * equilibrium in all age classes. Don't report any events. */
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


// ———  Set-up & tear-down  ———

Simulator::Simulator( util::Checksum ck, const scnXml::Scenario& scenario ) :
    simPeriodEnd(sim::zero()),
    totalSimDuration(sim::zero()),
    phase(STARTING_PHASE),
    workUnitIdentifier(0),
    cksum(ck)
{
    // ———  Initialise static data  ———
    
    const scnXml::Model& model = scenario.getModel();
    
    // 1) elements with no dependencies on other elements initialised here:
    sim::init( scenario );
    Parameters parameters( model.getParameters() );     // depends on nothing
    WithinHost::Genotypes::init( scenario );
    
    util::random::seed( model.getParameters().getIseed() );
    util::ModelOptions::init( model.getModelOptions() );
    
    // 2) elements depending on only elements initialised in (1):
    
    // Depends on parameters:
    if( scenario.getDiagnostics().present() ){
        WithinHost::diagnostics::init( parameters, scenario.getDiagnostics().get() );
    }
    
    // Survey init depends on diagnostics, monitoring:
    mon::initSurveyTimes( parameters, scenario, scenario.getMonitoring() );
    Population::init( parameters, scenario );
    
    // 3) elements depending on other elements; dependencies on (1) are not mentioned:
    
    // Transmission model initialisation depends on Transmission::PerHost and
    // genotypes (both from Human, from Population::init()) and
    // Monitoring::AgeGroup (from Surveys.init()):
    // Note: PerHost dependency can be postponed; it is only used to set adultAge
    population = auto_ptr<Population>(new Population( scenario.getEntomology(), scenario.getDemography().getPopSize() ));
    
    // Depends on transmission model (for species indexes):
    // MDA1D may depend on health system (too complex to verify)
    interventions::InterventionManager::init( scenario.getInterventions(), *population );
    
    // Depends on interventions, PK/PD (from population):
    Clinical::ClinicalModel::changeHS( scenario.getHealthSystem() );    // i.e. init health system
    
    // Depends on interventions:
    mon::initCohorts( scenario.getMonitoring() );
    
    // ———  End of static data initialisation  ———
    
    // Set work unit identifier, if we have one.
    if( scenario.getWuID().present() )
        workUnitIdentifier = scenario.getWuID().get();
    
    ifstream checkpointFile(CHECKPOINT,ios::in);
    // If not open, file doesn't exist (or is inaccessible)
    startedFromCheckpoint = checkpointFile.is_open();
}


// ———  run simulations  ———

void Simulator::start(const scnXml::Monitoring& monitoring){
    sim::time0 = sim::zero();
    sim::time1 = sim::zero();
    
    // Make sure warmup period is at least as long as a human lifespan, as the
    // length required by vector warmup, and is a whole number of years.
    SimTime humanWarmupLength = sim::maxHumanAge();
    if( humanWarmupLength < population->_transmissionModel->minPreinitDuration() ){
        cerr << "Warning: human life-span (" << humanWarmupLength.inYears();
        cerr << ") shorter than length of warm-up requested by" << endl;
        cerr << "transmission model ("
            << population->_transmissionModel->minPreinitDuration().inYears();
        cerr << "). Transmission may be unstable; perhaps use forced" << endl;
        cerr << "transmission (mode=\"forced\") or a longer life-span." << endl;
        humanWarmupLength = population->_transmissionModel->minPreinitDuration();
    }
    humanWarmupLength = sim::fromYearsI( ceil(humanWarmupLength.inYears()) );
    
    totalSimDuration = humanWarmupLength  // ONE_LIFE_SPAN
        + population->_transmissionModel->expectedInitDuration()
        // plus MAIN_PHASE: survey period plus one TS for last survey
        + mon::finalSurveyTime() + sim::oneTS();
    assert( totalSimDuration + sim::never() < sim::zero() );
    
    if (isCheckpoint()) {
        Continuous.init( monitoring, true );
        readCheckpoint();
    } else {
        Continuous.init( monitoring, false );
        population->createInitialHumans();
    }
    // Set to either a checkpointing time step or min int value. We only need to
    // set once, since we exit after a checkpoint triggered this way.
    SimTime testCheckpointTime = util::CommandLine::getNextCheckpointTime( sim::now() );
    SimTime testCheckpointDieTime = testCheckpointTime;        // kill program at same time
    
    // phase loop
    while (true){
        // loop for steps within a phase
        while (sim::now() < simPeriodEnd){
            // checkpoint
            if( util::BoincWrapper::timeToCheckpoint() || testCheckpointTime == sim::now() ){
                writeCheckpoint();
                util::BoincWrapper::checkpointCompleted();
            }
            if( testCheckpointDieTime == sim::now() ){
                throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
            }
            
            // Time step updates. Essentially, a time step is mid-day to
            // mid-day. sim::ts0() gives the date at the start of the step (in
            // internal units), and sim::ts1() the date at the end. Monitoring
            // and intervention deployment happen between updates, at time
            // sim::now().
            
            // Each step, monitoring (e.g. carrying out a survey) happens
            // first, reporting on the state at the start of the step (e.g.
            // patency diagnostics) or tallys of events which happened since
            // some point in the past (e.g. a previous survey).
            // Interventions are deployed next. For monitoring and intervention
            // deployment sim::ts0() and sim::ts1() should not be used.
            
            // Finally, Population::update1() runs, which updates both the
            // transmission model and the human population. During this time
            // sim::now() should not be used, however sim::ts0() equals its
            // last value while sim::ts1() is one greater.
            // Population::update1() also adds new humans, born at time
            //  sim::ts1(), to replace those lost.
            
            // do reporting (continuous and surveys)
            Continuous.update( *population );
            if( sim::intervNow() == mon::nextSurveyTime() ){
                population->newSurvey();
                mon::concludeSurvey();
            }
            
            // deploy interventions
            InterventionManager::deploy( *population );
            
            // update humans and mosquitoes
            
            // time1 is the time at the end of a time step, and mostly a
            // confusing relic, though sometimes useful.
            sim::time1 += sim::oneTS();
#ifndef NDEBUG
            sim::in_update = true;
#endif
            population->update1( humanWarmupLength );
#ifndef NDEBUG
            sim::in_update = false;
#endif
            sim::time0 += sim::oneTS();
            sim::interv_time += sim::oneTS();
            
            util::BoincWrapper::reportProgress(
                static_cast<double>(sim::now().raw()) /
                static_cast<double>(totalSimDuration.raw()) );
        }
        
        ++phase;        // advance to next phase
        if (phase == ONE_LIFE_SPAN) {
            simPeriodEnd = humanWarmupLength;
        } else if (phase == TRANSMISSION_INIT) {
            SimTime iterate = population->_transmissionModel->initIterate();
            if( iterate > sim::zero() ){
                simPeriodEnd += iterate;
                --phase;        // repeat phase
            } else {
                // nothing to do: start main phase immediately
            }
            // adjust estimation of final time step: end of current period + length of main phase
            totalSimDuration = simPeriodEnd + mon::finalSurveyTime() + sim::oneTS();
        } else if (phase == MAIN_PHASE) {
            // Start MAIN_PHASE:
            simPeriodEnd = totalSimDuration;
            sim::interv_time = sim::zero();
            population->preMainSimInit();
            population->newSurvey();       // Only to reset TransmissionModel::inoculationsPerAgeGroup
            mon::initMainSim();
        } else if (phase == END_SIM) {
            cerr << "sim end" << endl;
            break;
        }
        if (util::CommandLine::option (util::CommandLine::TEST_CHECKPOINTING)){
            // First of middle of next phase, or current value (from command line) triggers a checkpoint.
            SimTime phase_mid = sim::now() + (simPeriodEnd - sim::now()) * 0.5;
            // Don't checkpoint 0-length phases or do mid-phase checkpointing
            // when timed checkpoints were specified, and don't checkpoint
            // ONE_LIFE_SPAN phase if already past time humanWarmupLength:
            // these are extra transmission inits, and we don't want to
            // checkpoint every one of them.
            if( testCheckpointTime < sim::zero() && phase_mid > sim::now()
                && (phase != ONE_LIFE_SPAN || sim::now() < humanWarmupLength)
            ){
                testCheckpointTime = phase_mid;
                // Test checkpoint: die a bit later than checkpoint for better
                // resume testing (specifically, ctsout.txt).
                testCheckpointDieTime = testCheckpointTime + sim::oneTS() + sim::oneTS();
            }
        }
    }
    
    // Open a critical section; should prevent app kill while/after writing
    // output.txt, which we don't currently handle well.
    // Note: we don't end this critical section; we simply exit.
    util::BoincWrapper::beginCriticalSection();
    
    PopulationStats::print();
    
    population->flushReports();        // ensure all Human instances report past events
    mon::writeSurveyData();
    Continuous.finalise();
    
# ifdef OM_STREAM_VALIDATOR
    util::StreamValidator.saveStream();
# endif
}


// ———  checkpointing: set up read/write stream  ———

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

void Simulator::writeCheckpoint(){
    // We alternate between two checkpoints, in case program is closed while writing.
    const int NUM_CHECKPOINTS = 2;
    
    int oldCheckpointNum = 0, checkpointNum = 0;
    if (isCheckpoint()) {
        oldCheckpointNum = readCheckpointNum();
        // Get next checkpoint number:
        checkpointNum = mod_nn(oldCheckpointNum + 1, NUM_CHECKPOINTS);
    }
    
    {   // Open the next checkpoint file for writing:
        ostringstream name;
        name << CHECKPOINT << checkpointNum;
        //Writing checkpoint:
//      cerr << sim::now() << " WC: " << name.str();
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

void Simulator::readCheckpoint() {
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
  cerr << sim::now() << " RC" << endl;
  
  // On resume, write a checkpoint so we can tell whether we have identical checkpointed state
  if (util::CommandLine::option (util::CommandLine::TEST_DUPLICATE_CHECKPOINTS))
    writeCheckpoint();
}


// ———  checkpointing: Simulation data  ———

void Simulator::checkpoint (istream& stream, int checkpointNum) {
    try {
        util::checkpoint::header (stream);
        util::CommandLine::staticCheckpoint (stream);
        Population::staticCheckpoint (stream);
        Continuous & stream;
        mon::checkpoint( stream );
#       ifdef OM_STREAM_VALIDATOR
        util::StreamValidator & stream;
#       endif
        
        sim::interv_time & stream;
        simPeriodEnd & stream;
        totalSimDuration & stream;
        phase & stream;
        (*population) & stream;
        PopulationStats::staticCheckpoint( stream );
        InterventionManager::checkpoint( stream );
        InterventionManager::loadFromCheckpoint( *population, sim::interv_time );
        
        // read last, because other loads may use random numbers or expect time
        // to be negative
        sim::time0 & stream;
        sim::time1 & stream;
        util::random::checkpoint (stream, checkpointNum);
        
        // Check scenario.xml and checkpoint files correspond:
        int oldWUID = workUnitIdentifier;
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

void Simulator::checkpoint (ostream& stream, int checkpointNum) {
    util::checkpoint::header (stream);
    if (stream == NULL || !stream.good())
        throw util::checkpoint_error ("Unable to write to file");
    util::timer::startCheckpoint ();
    
    util::CommandLine::staticCheckpoint (stream);
    Population::staticCheckpoint (stream);
    Continuous & stream;
    mon::checkpoint( stream );
# ifdef OM_STREAM_VALIDATOR
    util::StreamValidator & stream;
# endif
    
    sim::interv_time & stream;
    simPeriodEnd & stream;
    totalSimDuration & stream;
    phase & stream;
    (*population) & stream;
    PopulationStats::staticCheckpoint( stream );
    InterventionManager::checkpoint( stream );
    
    sim::time0 & stream;
    sim::time1 & stream;
    util::random::checkpoint (stream, checkpointNum);
    workUnitIdentifier & stream;
    cksum & stream;
    
    util::timer::stopCheckpoint ();
    if (stream.fail())
        throw util::checkpoint_error ("stream write error");
}

}
