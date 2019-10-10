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
#include "mon/Continuous.h"
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
    using mon::Continuous;
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
    phase(STARTING_PHASE)
{
    // ———  Initialise static data  ———
    
    const scnXml::Model& model = scenario.getModel();
    
    // 1) elements with no dependencies on other elements initialised here:
    sim::init( scenario );  // also reads survey dates
    Parameters parameters( model.getParameters() );     // depends on nothing
    WithinHost::Genotypes::init( scenario );
    
    util::global_RNG.seed( model.getParameters().getIseed(), 721347520444481703 );
    // The master RNG is cryptographic with a hard-coded IV. Use of low
    // Hamming weight inputs (numbers close to 0) should not reduce quality.
    util::master_RNG.seed( model.getParameters().getIseed(), 0 );
    
    util::ModelOptions::init( model.getModelOptions() );
    
    // 2) elements depending on only elements initialised in (1):
    
    // Depends on parameters:
    WithinHost::diagnostics::init( parameters, scenario );
    
    // Reporting init depends on diagnostics, monitoring:
    mon::initReporting( scenario );
    Population::init( parameters, scenario );
    
    // 3) elements depending on other elements; dependencies on (1) are not mentioned:
    
    // Transmission model initialisation depends on Transmission::PerHost and
    // genotypes (both from Human, from Population::init()) and
    // mon::AgeGroup (from Surveys.init()):
    // Note: PerHost dependency can be postponed; it is only used to set adultAge
    population = unique_ptr<Population>(
            new Population( scenario.getDemography().getPopSize() ));
    uint64_t seed1 = util::master_RNG.gen_seed();
    uint64_t seed2 = util::master_RNG.gen_seed();
    transmission = unique_ptr<TransmissionModel>(
        TransmissionModel::createTransmissionModel(
            seed1, seed2,
            scenario.getEntomology(), population->size()) );
    
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
    sim::s_t0 = SimTime::zero();
    sim::s_t1 = SimTime::zero();
    
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
    
    m_estimatedEnd = humanWarmupLength  // ONE_LIFE_SPAN
        + transmission->expectedInitDuration()
        // plus MAIN_PHASE: survey period plus one TS for last survey
        + (sim::endDate() - sim::startDate())
        + SimTime::oneTS();
    assert( m_estimatedEnd + SimTime::never() < SimTime::zero() );
    
    if (isCheckpoint()) {
        Continuous.init( monitoring, true );
        readCheckpoint();
    } else {
        Continuous.init( monitoring, false );
        population->createInitialHumans();
        transmission->init2(*population);
    }
    
    int lastPercent = -1;	// last _integer_ percentage value
    
    // phase loop
    while (true){
        // loop for steps within a phase
        while (sim::now() < m_phaseEnd){
            int percent = (sim::now() * 100) / m_estimatedEnd;
            if( percent != lastPercent ){	// avoid huge amounts of output for performance/log-file size reasons
                lastPercent = percent;
                // \r cleans line. Then we print progress as a percentage.
                cerr << (boost::format("\r[%|3i|%%]\t") %percent) << flush;
            }
            
            // Monitoring. sim::now() gives time of end of last step,
            // and is when reporting happens in our time-series.
            Continuous.update( *population );
            if( sim::intervDate() == mon::nextSurveyDate() ){
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
            m_phaseEnd = humanWarmupLength;
            
        } else if (phase == TRANSMISSION_INIT) {
            // Start or continuation of transmission init cycle (after one life span)
            SimTime iterate = transmission->initIterate();
            if( iterate > SimTime::zero() ){
                m_phaseEnd += iterate;
                --phase;        // repeat phase
            } else {
                // nothing to do: start main phase immediately
            }
            // adjust estimation of final time step: end of current period + length of main phase
            m_estimatedEnd = m_phaseEnd
                + (sim::endDate() - sim::startDate())
                + SimTime::oneTS();
            
        } else if (phase == MAIN_PHASE) {
            // Start MAIN_PHASE:
            m_phaseEnd = m_estimatedEnd;
            sim::s_interv = SimTime::zero();
            population->preMainSimInit();
            transmission->summarize();    // Only to reset TransmissionModel::inoculationsPerAgeGroup
            mon::initMainSim();
            
        } else if (phase == END_SIM) {
            cerr << "sim end" << endl;
            break;
        }
        
        if (phase == MAIN_PHASE && util::CommandLine::option (util::CommandLine::CHECKPOINT)){
            writeCheckpoint();
            if( util::CommandLine::option (util::CommandLine::CHECKPOINT_STOP) ){
                throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
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
        name << CHECKPOINT << checkpointNum << ".gz";
        //Writing checkpoint:
        ogzstream out(name.str().c_str(), ios::out | ios::binary);
        checkpoint (out);
        out.close();
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
    if( oldCheckpointNum != checkpointNum ){
        ostringstream name;
        name << CHECKPOINT << oldCheckpointNum << ".gz";
        ofstream out(name.str().c_str(), ios::out | ios::binary);
        out.close();
    }
//     cerr << " OK" << endl;
}

void Simulator::readCheckpoint() {
    int checkpointNum = readCheckpointNum();
    
    // Open the latest file
    ostringstream name;
    name << CHECKPOINT << checkpointNum << ".gz";
    igzstream in(name.str().c_str(), ios::in | ios::binary);
    //Note: gzstreams are considered "good" when file not open!
    if ( !( in.good() && in.rdbuf()->is_open() ) )
        throw util::checkpoint_error ("Unable to read file");
    checkpoint (in);
    in.close();
  
    cerr << sim::now().inSteps() << "t loaded checkpoint" << endl;
}


// ———  checkpointing: Simulation data  ———

void Simulator::checkpoint (istream& stream) {
    try {
        util::checkpoint::header (stream);
        util::CommandLine::staticCheckpoint (stream);
        Population::staticCheckpoint (stream);
        Continuous & stream;
        mon::checkpoint( stream );
#       ifdef OM_STREAM_VALIDATOR
        util::StreamValidator & stream;
#       endif
        
        sim::s_interv & stream;
        m_phaseEnd & stream;
        m_estimatedEnd & stream;
        phase & stream;
        transmission & stream;
        population->checkpoint(stream);
        InterventionManager::checkpoint( stream );
        InterventionManager::loadFromCheckpoint( *population, *transmission );
        
        // read last, because other loads may use random numbers or expect time
        // to be negative
        sim::s_t0 & stream;
        sim::s_t1 & stream;
        util::global_RNG.checkpoint(stream);
        util::master_RNG.checkpoint(stream);
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

void Simulator::checkpoint (ostream& stream) {
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
    
    sim::s_interv & stream;
    m_phaseEnd & stream;
    m_estimatedEnd & stream;
    phase & stream;
    transmission & stream;
    population->checkpoint(stream);
    InterventionManager::checkpoint( stream );
    
    sim::s_t0 & stream;
    sim::s_t1 & stream;
    util::global_RNG.checkpoint (stream);
    util::master_RNG.checkpoint(stream);
    
    util::timer::stopCheckpoint ();
    if (stream.fail())
        throw util::checkpoint_error ("stream write error");
}

}
