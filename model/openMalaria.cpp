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

#include "util/DocumentLoader.h"

#include "Global.h"
#include "Transmission/transmission.h"
#include "Parameters.h"
#include "Clinical/ClinicalModel.h"
#include "mon/Continuous.h"
#include "interventions/InterventionManager.hpp"
#include "Population.h"
#include "WithinHost/WHInterface.h"
#include "WithinHost/Diagnostic.h"
#include "WithinHost/Genotypes.h"
#include "mon/management.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "schema/scenario.h"

#include <fstream>
#include <gzstream/gzstream.h>

#include <cstdio>
#include <cerrno>

namespace scnXml{
    class Monitoring;
    class Scenario;
}

namespace OM{
    using mon::Continuous;
    using interventions::InterventionManager;
    using Transmission::TransmissionModel;
}

using namespace OM;

/** @brief checkpointing functions
*
* readCheckpoint/writeCheckpoint prepare to read/write the file,
* and read/write read and write the actual data. */
//@{
int readCheckpointNum (const string &checkpointFileName)
{
    ifstream checkpointFile;
    checkpointFile.open(checkpointFileName, fstream::in);
    int checkpointNum=0;
    checkpointFile >> checkpointNum;
    checkpointFile.close();
    if (!checkpointFile)
        throw util::checkpoint_error ("error reading from file \"checkpoint\"");
    return checkpointNum;
}

void checkpoint (istream& stream, SimTime &endTime, SimTime &estEndTime, Population &population, TransmissionModel &transmission)
{
    try
    {
        util::checkpoint::header (stream);
        util::CommandLine::staticCheckpoint (stream);
        Population::staticCheckpoint (stream);
        Continuous & stream;
        mon::checkpoint( stream );
#       ifdef OM_STREAM_VALIDATOR
        util::StreamValidator & stream;
#       endif
        
        sim::s_interv & stream;
        endTime & stream;
        estEndTime & stream;
        transmission & stream;
        population.checkpoint(stream);
        InterventionManager::checkpoint(stream);
        InterventionManager::loadFromCheckpoint(population, transmission);
        
        // read last, because other loads may use random numbers or expect time
        // to be negative
        sim::s_t0 & stream;
        sim::s_t1 & stream;
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

void checkpoint (ostream& stream, SimTime &endTime, SimTime &estEndTime, Population &population, TransmissionModel &transmission) {
    util::checkpoint::header (stream);
    if (!stream.good())
        throw util::checkpoint_error ("Unable to write to file");

    util::CommandLine::staticCheckpoint (stream);
    Population::staticCheckpoint (stream);
    Continuous & stream;
    mon::checkpoint( stream );
# ifdef OM_STREAM_VALIDATOR
    util::StreamValidator & stream;
# endif
    
    sim::s_interv & stream;
    endTime & stream;
    estEndTime & stream;
    transmission & stream;
    population.checkpoint(stream);
    InterventionManager::checkpoint( stream );
    
    sim::s_t0 & stream;
    sim::s_t1 & stream;
    util::master_RNG.checkpoint(stream);
    
    if (stream.fail())
        throw util::checkpoint_error ("stream write error");
}

void writeCheckpoint(const bool startedFromCheckpoint, const string &checkpointFileName, SimTime &endTime, SimTime &estEndTime, Population &population, TransmissionModel &transmission)
{
    // We alternate between two checkpoints, in case program is closed while writing.
    const int NUM_CHECKPOINTS = 2;
    
    int oldCheckpointNum = 0, checkpointNum = 0;
    if (startedFromCheckpoint)
    {
        oldCheckpointNum = readCheckpointNum(checkpointFileName);
        checkpointNum = mod_nn(oldCheckpointNum + 1, NUM_CHECKPOINTS); // Get next checkpoint number:
    }
    
    {   // Open the next checkpoint file for writing:
        ostringstream name;
        name << checkpointFileName << checkpointNum << ".gz";
        ogzstream out(name.str().c_str(), ios::out | ios::binary);
        checkpoint (out, endTime, estEndTime, population, transmission);
        out.close();
    }
    
    {   // Indicate which is the latest checkpoint file.
        ofstream checkpointFile;
        checkpointFile.open(checkpointFileName,ios::out);
        checkpointFile << checkpointNum;
        checkpointFile.close();
        if (!checkpointFile)
            throw util::checkpoint_error ("error writing to file \"checkpoint\"");
    }

    // Truncate the old checkpoint to save disk space, when it existed
    if( oldCheckpointNum != checkpointNum ){
        ostringstream name;
        name << checkpointFileName << oldCheckpointNum << ".gz";
        ofstream out(name.str().c_str(), ios::out | ios::binary);
        out.close();
    }
}

void readCheckpoint(const string &checkpointFileName, SimTime &endTime, SimTime &estEndTime, Population &population, TransmissionModel &transmission)
{
    int checkpointNum = readCheckpointNum(checkpointFileName);
    
    // Open the latest file
    ostringstream name;
    name << checkpointFileName << checkpointNum << ".gz";
    igzstream in(name.str().c_str(), ios::in | ios::binary);
    //Note: gzstreams are considered "good" when file not open!
    if ( !( in.good() && in.rdbuf()->is_open() ) )
        throw util::checkpoint_error ("Unable to read file");
    checkpoint (in, endTime, estEndTime, population, transmission);
    in.close();
  
    cerr << sim::inSteps(sim::now()) << "t loaded checkpoint" << endl;
}

// Internal simulation loop
void loop(const SimTime humanWarmupLength, Population &population, TransmissionModel &transmission, SimTime &endTime, SimTime &estEndTime, int lastPercent)
{
    while (sim::now() < endTime)
    {        
        // Monitoring. sim::now() gives time of end of last step,
        // and is when reporting happens in our time-series.
        Continuous.update( population );
        if( sim::intervDate() == mon::nextSurveyDate() ){
            population.newSurvey();
            transmission.summarize();
            mon::concludeSurvey();
        }
        
        // Deploy interventions, at time sim::now().
        InterventionManager::deploy( population, transmission );
        
        // Time step updates. Time steps are mid-day to mid-day.
        // sim::ts0() gives the date at the start of the step, sim::ts1() the date at the end.
        sim::start_update();
        
        // This should be called before humans contract new infections in the simulation step.
        // This needs the whole population (it is an approximation before all humans are updated).
        transmission.vectorUpdate (population);
        
        population.update(transmission, humanWarmupLength);
        
        // Doesn't matter whether non-updated humans are included (value isn't used
        // before all humans are updated).
        transmission.update(population);
        
        sim::end_update();

        int percent = (sim::now() * 100) / estEndTime;
        if( percent != lastPercent ){   // avoid huge amounts of output for performance/log-file size reasons
            lastPercent = percent;
            // \r cleans line. Then we print progress as a percentage.
            cerr << "\r" << percent << "%\t" << flush;
        }

        if( errno != 0 )
        {
           char err[256];
           sprintf(err, "t = %d Please report! Error: ", sim::now());
           std::perror(err);
           errno = 0;
        }
    }
}

/// main() — loads scenario XML and runs simulation
int main(int argc, char* argv[])
{
    int exitStatus = EXIT_SUCCESS;
    string scenarioFile;
    bool startedFromCheckpoint;
    string checkpointFileName;
    SimTime estEndTime, endTime;
    
    try {
        util::set_gsl_handler();        // init
        
        scenarioFile = util::CommandLine::parse (argc, argv);   // parse arguments
        
        // Load the scenario document:
        scenarioFile = util::CommandLine::lookupResource (scenarioFile);
        util::DocumentLoader documentLoader;
        documentLoader.loadDocument(scenarioFile);
        
        const scnXml::Scenario &scenario = documentLoader.document();
        const scnXml::Model &model = scenario.getModel();
        const scnXml::Monitoring &monitoring = documentLoader.document().getMonitoring();
        
        // 1) elements with no dependencies on other elements initialised here:
        sim::init( scenario );  // also reads survey dates
        Parameters parameters( model.getParameters() );     // depends on nothing
        WithinHost::Genotypes::init( scenario );
        
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
        size_t popSize = scenario.getDemography().getPopSize();
        std::unique_ptr<Population> population = unique_ptr<Population>(new Population( popSize ));
        std::unique_ptr<TransmissionModel> transmission = unique_ptr<TransmissionModel>(Transmission::createTransmissionModel(scenario.getEntomology(), popSize));
        
        // Depends on transmission model (for species indexes):
        // MDA1D may depend on health system (too complex to verify)
        interventions::InterventionManager::init(scenario.getInterventions(), *population, *transmission );
        
        // Depends on interventions, PK/PD (from humanPop):
        Clinical::ClinicalModel::setHS( scenario.getHealthSystem() );
        
        // Depends on interventions:
        mon::initCohorts( scenario.getMonitoring() );
        
        // ———  End of static data initialisation  ———
        checkpointFileName = util::CommandLine::getCheckpointName();

        if(checkpointFileName == "")
            checkpointFileName = "checkpoint";

        if(util::CommandLine::option(util::CommandLine::CHECKPOINT))
        {
            ifstream checkpointFile(checkpointFileName,ios::in);
            // If not open, file doesn't exist (or is inaccessible)
            startedFromCheckpoint = checkpointFile.is_open();

            // Cleanup errno in file doesn't exist
            if(startedFromCheckpoint == false)
                errno = 0;
        }
        else
            startedFromCheckpoint = false;

        sim::s_t0 = sim::zero();
        sim::s_t1 = sim::zero();
        
        // Make sure warmup period is at least as long as a human lifespan, as the
        // length required by vector warmup, and is a whole number of years.
        SimTime humanWarmupLength = sim::maxHumanAge();
        if( humanWarmupLength < transmission->minPreinitDuration() ){
            cerr << "Warning: human life-span (" << sim::inYears(humanWarmupLength);
            cerr << ") shorter than length of warm-up requested by" << endl;
            cerr << "transmission model ("
                << sim::inYears(transmission->minPreinitDuration());
            cerr << "). Transmission may be unstable; perhaps use forced" << endl;
            cerr << "transmission (mode=\"forced\") or a longer life-span." << endl;
            humanWarmupLength = transmission->minPreinitDuration();
        }
        humanWarmupLength = sim::fromYearsI( static_cast<int>(ceil(sim::inYears(humanWarmupLength))) );
        
        estEndTime = humanWarmupLength  // ONE_LIFE_SPAN
            + transmission->expectedInitDuration()
            // plus MAIN_PHASE: survey period plus one TS for last survey
            + (sim::endDate() - sim::startDate())
            + sim::oneTS();
        assert( estEndTime + sim::never() < sim::zero() );
        
        bool skipWarmup = false;
        if (startedFromCheckpoint)
        {
            Continuous.init( monitoring, true );
            readCheckpoint(checkpointFileName, endTime, estEndTime, *population, *transmission);
            skipWarmup = true;
        }
        else
        {
            Continuous.init( monitoring, false );
            population->createInitialHumans();
            transmission->init2(*population);
        }
        
        int lastPercent = -1;   // last _integer_ percentage value
        
        if(!skipWarmup)
        {
            /** Warm-up phase: 
             * Run the simulation using the equilibrium inoculation rates over one
             * complete lifespan (sim::maxHumanAge()) to reach immunological
             * equilibrium in all age classes. Don't report any events. */
            endTime = humanWarmupLength;
            loop(humanWarmupLength, *population, *transmission, endTime, estEndTime, lastPercent);

            // Transmission init phase
            SimTime iterate = transmission->initIterate();
            while(iterate > sim::zero())
            {
                endTime += iterate;
                // adjust estimation of final time step: end of current period + length of main phase
                estEndTime = endTime + (sim::endDate() - sim::startDate()) + sim::oneTS();
                loop(humanWarmupLength, *population, *transmission, endTime, estEndTime, lastPercent);
                iterate = transmission->initIterate();
            }

            // Main phase
            //! This procedure starts with the current state of the simulation 
            /*! It continues updating assuming:
                (i)         the default (exponential) demographic model
                (ii)        the entomological input defined by the EIRs in intEIR()
                (iii)       the intervention packages defined in Intervention()
                (iv)        the survey times defined in Survey() */
            endTime = estEndTime;
            sim::s_interv = sim::zero();
            population->preMainSimInit();
            transmission->summarize();    // Only to reset TransmissionModel::inoculationsPerAgeGroup
            mon::initMainSim();

            if(util::CommandLine::option (util::CommandLine::CHECKPOINT))
            {
                writeCheckpoint(startedFromCheckpoint, checkpointFileName, endTime, estEndTime, *population, *transmission);
                if( util::CommandLine::option (util::CommandLine::CHECKPOINT_STOP) )
                    throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
            }
        }

        loop(humanWarmupLength, *population, *transmission, endTime, estEndTime, lastPercent);
       
        cerr << '\r' << flush;  // clean last line of progress-output
        
        population->flushReports();        // ensure all Human instances report past events
        mon::writeSurveyData();
        
    # ifdef OM_STREAM_VALIDATOR
        util::StreamValidator.saveStream();
    # endif
        
        // simulation's destructor runs
    } catch (const OM::util::cmd_exception& e) {
        if( e.getCode() == 0 ){
            // this is not an error, but exiting due to command line
            cerr << e.what() << "; exiting..." << endl;
        }else{
            cerr << "Command-line error: "<<e.what();
            exitStatus = e.getCode();
        }
    } catch (const ::xsd::cxx::tree::exception<char>& e) {
        cerr << "XSD error: " << e.what() << '\n' << e << endl;
        exitStatus = OM::util::Error::XSD;
    } catch (const OM::util::checkpoint_error& e) {
        cerr << "Checkpoint error: " << e.what() << endl;
        cerr << e << flush;
        exitStatus = e.getCode();
    } catch (const OM::util::traced_exception& e) {
        cerr << "Code error: " << e.what() << endl;
        cerr << e << flush;
        cerr << "This is likely an error in the C++ code. Please report!" << endl;
        exitStatus = e.getCode();
    } catch (const OM::util::xml_scenario_error& e) {
        cerr << "Error: " << e.what() << endl;
        cerr << "In: " << scenarioFile << endl;
        exitStatus = e.getCode();
    } catch (const OM::util::base_exception& e) {
        cerr << "Error: " << e.message() << endl;
        exitStatus = e.getCode();
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        exitStatus = EXIT_FAILURE;
    } catch (...) {
        cerr << "Unknown error" << endl;
        exitStatus = EXIT_FAILURE;
    }
    
    // If we get to here, we already know an error occurred.
    if( errno != 0 )
        std::perror( "OpenMalaria" );
    
    return exitStatus;
}
