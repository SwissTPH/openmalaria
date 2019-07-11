/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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
#include "Parameters.h"
#include "Clinical/ClinicalModel.h"
#include "Monitoring/Continuous.h"
#include "interventions/InterventionManager.hpp"
#include "Population.h"
#include "WithinHost/WHInterface.h"
#include "WithinHost/Diagnostic.h"
#include "WithinHost/Genotypes.h"
#include "mon/management.h"
#include "util/timer.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "schema/scenario.h"

#include <fstream>
#include <gzstream/gzstream.h>
#include <boost/format.hpp>


namespace OM {
    using Monitoring::Continuous;
    using interventions::InterventionManager;
    using Transmission::TransmissionModel;

bool Simulator::startedFromCheckpoint;  // static

const char* CHECKPOINT = "checkpoint";

std::unique_ptr<Population> population;
std::unique_ptr<TransmissionModel> transmission;

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

Simulator::Simulator( const scnXml::Scenario& scenario ) :
    simPeriodEnd(SimTime::zero()),
    totalSimDuration(SimTime::zero()),
    phase(STARTING_PHASE)
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
    WithinHost::diagnostics::init( parameters, scenario );
    
    // Survey init depends on diagnostics, monitoring:
    mon::initSurveyTimes( parameters, scenario, scenario.getMonitoring() );
    Population::init( parameters, scenario );
    
    // 3) elements depending on other elements; dependencies on (1) are not mentioned:
    
    // Transmission model initialisation depends on Transmission::PerHost and
    // genotypes (both from Human, from Population::init()) and
    // Monitoring::AgeGroup (from Surveys.init()):
    // Note: PerHost dependency can be postponed; it is only used to set adultAge
    population = unique_ptr<Population>(
            new Population( scenario.getDemography().getPopSize() ));
    transmission = unique_ptr<TransmissionModel>(
            TransmissionModel::createTransmissionModel(scenario.getEntomology(), population->size()) );
    
    // Depends on transmission model (for species indexes):
    // MDA1D may depend on health system (too complex to verify)
    interventions::InterventionManager::init( scenario.getInterventions(), *transmission );
    
    // Depends on interventions, PK/PD (from humanPop):
    Clinical::ClinicalModel::setHS( scenario.getHealthSystem() );
    
    // Depends on interventions:
    mon::initCohorts( scenario.getMonitoring() );
    
    // ———  End of static data initialisation  ———
    
    ifstream checkpointFile(CHECKPOINT,ios::in);
    // If not open, file doesn't exist (or is inaccessible)
    startedFromCheckpoint = checkpointFile.is_open();
}


// ———  run simulations  ———

void Simulator::start(const scnXml::Monitoring& monitoring){
    sim::time0 = SimTime::zero();
    sim::time1 = SimTime::zero();
    
    // Make sure warmup period is at least as long as a human lifespan, as the
    // length required by vector warmup, and is a whole number of years.
    SimTime humanWarmupLength = sim::maxHumanAge();
    if( humanWarmupLength < transmission->minPreinitDuration() ){
        cerr << "Warning: human life-span (" << humanWarmupLength.inYears();
        cerr << ") shorter than length of warm-up requested by" << endl;
        cerr << "transmission model ("
            << transmission->minPreinitDuration().inYears();
        cerr << "). Transmission may be unstable; perhaps use forced" << endl;
        cerr << "transmission (mode=\"forced\") or a longer life-span." << endl;
        humanWarmupLength = transmission->minPreinitDuration();
    }
    humanWarmupLength = SimTime::fromYearsI( static_cast<int>(ceil(humanWarmupLength.inYears())) );
    
    totalSimDuration = humanWarmupLength  // ONE_LIFE_SPAN
        + transmission->expectedInitDuration()
        // plus MAIN_PHASE: survey period plus one TS for last survey
        + mon::finalSurveyTime() + SimTime::oneTS();
    assert( totalSimDuration + SimTime::never() < SimTime::zero() );
    
    if (isCheckpoint()) {
        Continuous.init( monitoring, true );
        readCheckpoint();
    } else {
        Continuous.init( monitoring, false );
        population->createInitialHumans();
        transmission->init2(*population);
    }
    // Set to either a checkpointing time step or min int value. We only need to
    // set once, since we exit after a checkpoint triggered this way.
    SimTime testCheckpointTime = util::CommandLine::getNextCheckpointTime( sim::now() );
    SimTime testCheckpointDieTime = testCheckpointTime;        // kill program at same time
    
    int lastPercent = -1;	// last _integer_ percentage value
    
    // phase loop
    while (true){
        // loop for steps within a phase
        while (sim::now() < simPeriodEnd){
            int percent = (sim::now() * 100) / totalSimDuration;
            if( percent != lastPercent ){	// avoid huge amounts of output for performance/log-file size reasons
                lastPercent = percent;
                // \r cleans line. Then we print progress as a percentage.
                cerr << (boost::format("\r[%|3i|%%]\t") %percent) << flush;
            }
            
            if( testCheckpointTime == sim::now() ){
                writeCheckpoint();
            }
            if( testCheckpointDieTime == sim::now() ){
                throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
            }
            
            // Monitoring. sim::now() gives time of end of last step,
            // and is when reporting happens in our time-series.
            Continuous.update( *population );
            if( sim::intervNow() == mon::nextSurveyTime() ){
                population->newSurvey();
                transmission->summarize();
                mon::concludeSurvey();
            }
            
            // Deploy interventions, at time sim::now().
            InterventionManager::deploy( *population, *transmission );
            
            // Time step updates. Time steps are mid-day to mid-day.
            // sim::ts0() gives the date at the start of the step, sim::ts1() the date at the end.
            sim::start_update();
            
            // This should be called before humans contract new infections in the simulation step.
            // This needs the whole population (it is an approximation before all humans are updated).
            transmission->vectorUpdate (*population);
            
            population->update(*transmission, humanWarmupLength);
            
            // Doesn't matter whether non-updated humans are included (value isn't used
            // before all humans are updated).
            transmission->update(*population);
            
            sim::end_update();
        }
        
        ++phase;        // advance to next phase
        if (phase == ONE_LIFE_SPAN) {
            // Start human warm-up
            simPeriodEnd = humanWarmupLength;
            
        } else if (phase == TRANSMISSION_INIT) {
            // Start or continuation of transmission init cycle (after one life span)
            SimTime iterate = transmission->initIterate();
            if( iterate > SimTime::zero() ){
                simPeriodEnd += iterate;
                --phase;        // repeat phase
            } else {
                // nothing to do: start main phase immediately
            }
            // adjust estimation of final time step: end of current period + length of main phase
            totalSimDuration = simPeriodEnd + mon::finalSurveyTime() + SimTime::oneTS();
            
        } else if (phase == MAIN_PHASE) {
            // Start MAIN_PHASE:
            simPeriodEnd = totalSimDuration;
            sim::interv_time = SimTime::zero();
            population->preMainSimInit();
            transmission->summarize();    // Only to reset TransmissionModel::inoculationsPerAgeGroup
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
            if( testCheckpointTime < SimTime::zero() && phase_mid > sim::now()
                && (phase != ONE_LIFE_SPAN || sim::now() < humanWarmupLength)
            ){
                testCheckpointTime = phase_mid;
                // Test checkpoint: die a bit later than checkpoint for better
                // resume testing (specifically, ctsout.txt).
                testCheckpointDieTime = testCheckpointTime + SimTime::oneTS() + SimTime::oneTS();
            }
        }
    }
    
    cerr << '\r' << flush;	// clean last line of progress-output
    
    population->flushReports();        // ensure all Human instances report past events
    mon::writeSurveyData();
    
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
  cerr << sim::now().inSteps() << "t RC" << endl;
  
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
        transmission & stream;
        population->checkpoint(stream);
        InterventionManager::checkpoint( stream );
        InterventionManager::loadFromCheckpoint( *population, *transmission, sim::interv_time );
        
        // read last, because other loads may use random numbers or expect time
        // to be negative
        sim::time0 & stream;
        sim::time1 & stream;
        util::random::checkpoint (stream, checkpointNum);
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
    if (!stream.good())
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
    transmission & stream;
    population->checkpoint(stream);
    InterventionManager::checkpoint( stream );
    
    sim::time0 & stream;
    sim::time1 & stream;
    util::random::checkpoint (stream, checkpointNum);
    
    util::timer::stopCheckpoint ();
    if (stream.fail())
        throw util::checkpoint_error ("stream write error");
}

}
