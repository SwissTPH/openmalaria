/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
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

#include "Global.h"

#include "Transmission/transmission.h"
#include "Clinical/ClinicalModel.h"
#include "interventions/InterventionManager.hpp"

#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include "util/DocumentLoader.h"

#include "mon/Continuous.h"
#include "mon/management.h"

#include "schema/scenario.h"

#include "Population.h"
#include "Parameters.h"
#include "checkpoint.h"

#include "Host/NeonatalMortality.h"

#include <cerrno>

namespace OM {
    using mon::Continuous;
    using interventions::InterventionManager;
    using Transmission::TransmissionModel;
}

using namespace OM;

void print_progress(int lastPercent, SimTime &estEndTime)
{
    int percent = (sim::now() * 100) / estEndTime;
    if( percent != lastPercent ) {   // avoid huge amounts of output for performance/log-file size reasons
        lastPercent = percent;
        cerr << "\r" << percent << "%\t" << flush;
    }
}

void print_errno()
{
    if( errno != 0 )
    {
       char err[256];
       sprintf(err, "t = %d Please report! Error: ", int(sim::now()));
       std::perror(err);
       errno = 0;
    }
}

// Internal simulation loop
void loop(const SimTime humanWarmupLength, Population &population, TransmissionModel &transmission, SimTime &endTime, SimTime &estEndTime, int lastPercent)
{
    while (sim::now() < endTime)
    {
        if (util::CommandLine::option(util::CommandLine::VERBOSE) && sim::intervDate() > 0)
            cout << "Time step: " << sim::now() / sim::oneTS() << ", internal days: " << sim::now() << " | " << estEndTime << ", Intervention Date: " << sim::intervDate() << endl;

        // Monitoring. sim::now() gives time of end of last step,
        // and is when reporting happens in our time-series.
        Continuous.update( population );
        if( sim::intervDate() == mon::nextSurveyDate() ){
            for(Host::Human &human : population.humans)
                human.summarize();
            transmission.summarize();
            mon::concludeSurvey();
        }
        
        // Deploy interventions, at time sim::now().
        InterventionManager::deploy( population.humans, transmission );
        
        // Time step updates. Time steps are mid-day to mid-day.
        // sim::ts0() gives the date at the start of the step, sim::ts1() the date at the end.
        sim::start_update();

        // This should be called before humans contract new infections in the simulation step.
        // This needs the whole population (it is an approximation before all humans are updated).
        transmission.vectorUpdate(population.humans);
        
        // NOTE: no neonatal mortalities will occur in the first 20 years of warmup
        // (until humans old enough to be pregnate get updated and can be infected).
        Host::NeonatalMortality::update (population.humans);
        
        for (Host::Human& human : population.humans)
            if (human.getDateOfBirth() + sim::maxHumanAge() >= humanWarmupLength) // this is last time of possible update
                human.update(transmission);
       
        population.regularize();
        
        // Doesn't matter whether non-updated humans are included (value isn't used
        // before all humans are updated).
        transmission.updateKappa(population.humans);

        sim::end_update();

        print_progress(lastPercent, estEndTime);
        print_errno();
    }
}

/// main() — loads scenario XML and runs simulation
int main(int argc, char* argv[])
{
    int exitStatus = EXIT_SUCCESS;
    bool startedFromCheckpoint;

    string scenarioFile;
    string checkpointFileName;
    SimTime estEndTime, endTime;
    
    try {
        util::set_gsl_handler();
        
        scenarioFile = util::CommandLine::parse (argc, argv);
        unique_ptr<scnXml::Scenario> scenario = util::loadScenario(scenarioFile);

        const scnXml::Model &model = scenario->getModel();
        const scnXml::Monitoring &monitoring = scenario->getMonitoring();
        
        // 1) elements with no dependencies on other elements initialised here:
        sim::init( *scenario );  // also reads survey dates
        Parameters parameters( model.getParameters() );     // depends on nothing
        WithinHost::Genotypes::init( *scenario );
        
        util::master_RNG.seed( model.getParameters().getIseed(), 0 ); // Init RNG with Iseed
        util::ModelOptions::init( model.getModelOptions() );
        
        // 2) elements depending on only elements initialised in (1):
        WithinHost::diagnostics::init( parameters, *scenario ); // Depends on Parameters
        mon::initReporting( *scenario ); // Reporting init depends on diagnostics and monitoring
        Host::Human::init( parameters, *scenario );
        Host::NeonatalMortality::init( scenario->getModel().getClinical() );
        AgeStructure::init( scenario->getDemography() );
        
        // 3) elements depending on other elements; dependencies on (1) are not mentioned:
        // Transmission model initialisation depends on Transmission::PerHost and
        // genotypes (both from Human, from Population::init()) and
        // mon::AgeGroup (from Surveys.init()):
        // Note: PerHost dependency can be postponed; it is only used to set adultAge
        size_t popSize = scenario->getDemography().getPopSize();
        std::unique_ptr<Population> population = unique_ptr<Population>(population::createPopulation(popSize));
        std::unique_ptr<TransmissionModel> transmission = unique_ptr<TransmissionModel>(Transmission::createTransmissionModel(scenario->getEntomology(), popSize));

        // Depends on transmission model (for species indexes):
        // MDA1D may depend on health system (too complex to verify)
        interventions::InterventionManager::init(scenario->getInterventions(), *population, *transmission );
        Clinical::ClinicalModel::setHS( scenario->getHealthSystem() ); // Depends on interventions, PK/PD (from humanPop)
        mon::initCohorts( scenario->getMonitoring() ); // Depends on interventions
        
        // ———  End of static data initialisation  ———
        checkpointFileName = util::CommandLine::getCheckpointName();

        if(checkpointFileName == "")
            checkpointFileName = "checkpoint";

        if(util::CommandLine::option(util::CommandLine::CHECKPOINT))
        {
            ifstream checkpointFile(checkpointFileName,ios::in);
            startedFromCheckpoint = checkpointFile.is_open();
            if(startedFromCheckpoint == false)
                errno = 0; // Cleanup errno if file doesn't exist
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
            cerr << "transmission model (" << sim::inYears(transmission->minPreinitDuration());
            cerr << "). Transmission may be unstable; perhaps use forced" << endl;
            cerr << "transmission (mode=\"forced\") or a longer life-span." << endl;
            humanWarmupLength = transmission->minPreinitDuration();
        }
        humanWarmupLength = sim::fromYearsI( static_cast<int>(ceil(sim::inYears(humanWarmupLength))) );
        
        estEndTime = humanWarmupLength + transmission->expectedInitDuration() + (sim::endDate() - sim::startDate()) + sim::oneTS();
        assert( estEndTime + sim::never() < sim::zero() );
        
        int lastPercent = -1; // last _integer_ percentage value
        
        if (startedFromCheckpoint)
        {
            Continuous.init( monitoring, true );
            readCheckpoint(checkpointFileName, endTime, estEndTime, *population, *transmission);
        }
        else
        {
            Continuous.init( monitoring, false );
            population->createInitialHumans();
            transmission->init2(population->humans);
        
            // ofstream humanFile ("humans.txt");
            // if (humanFile.is_open())
            // {
            //     for(size_t i=0; i<population->getHumans()[0].perHostTransmission.speciesData.size(); i++)
            //     {
            //         for(const Host::Human &human : population->getHumans())
            //         {
            //             humanFile << human.perHostTransmission.speciesData[i].getEntoAvailability() << " ";
            //         }
            //         humanFile << endl;
            //     }
            //     humanFile.close();
            // }

            /** Warm-up phase: 
             * Run the simulation using the equilibrium inoculation rates over one
             * complete lifespan (sim::maxHumanAge()) to reach immunological
             * equilibrium in all age classes. Don't report any events. */
            endTime = humanWarmupLength;
            if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Starting Warmup..." << endl;
            loop(humanWarmupLength, *population, *transmission, endTime, estEndTime, lastPercent);
            if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Finishing Warmup..." << endl;


            /** Transmission init phase:
             * Fit the emergence rate to the input EIR */
            SimTime iterate = transmission->initIterate();
            while(iterate > sim::zero())
            {
                endTime = endTime + iterate;
                // adjust estimation of final time step: end of current period + length of main phase
                estEndTime = endTime + (sim::endDate() - sim::startDate()) + sim::oneTS();
                if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Starting EIR Calibration..." << endl;
                loop(humanWarmupLength, *population, *transmission, endTime, estEndTime, lastPercent);
                if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Finishing EIR Calibration..." << endl;
                iterate = transmission->initIterate();
            }

            /** Main phase:
             * This procedure starts with the current state of the simulation 
             * It continues updating assuming:
             * (i)         the default (exponential) demographic model
             * (ii)        the entomological input defined by the EIRs in intEIR()
             * (iii)       the intervention packages defined in Intervention()
             * (iv)        the survey times defined in Survey() */
            // reset endTime and estEndTime to their exact value after initIterate()
            estEndTime = endTime = endTime + (sim::endDate() - sim::startDate()) + sim::oneTS();
            sim::s_interv = sim::zero();
            Host::InfectionIncidenceModel::preMainSimInit();
            Clinical::InfantMortality::preMainSimInit();
            WithinHost::Genotypes::preMainSimInit();
            population->recentBirths = 0;
            transmission->summarize(); // Only to reset TransmissionModel::inoculationsPerAgeGroup
            mon::initMainSim();

            if(util::CommandLine::option (util::CommandLine::CHECKPOINT))
            {
                writeCheckpoint(startedFromCheckpoint, checkpointFileName, endTime, estEndTime, *population, *transmission);
                if( util::CommandLine::option (util::CommandLine::CHECKPOINT_STOP) )
                    throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
            }
        }

        // Main phase loop
        if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Starting Intervention period..." << endl;
        loop(humanWarmupLength, *population, *transmission, endTime, estEndTime, lastPercent);
        if (util::CommandLine::option(util::CommandLine::VERBOSE)) cout << "Finishing Intervention period..." << endl;
       
        cerr << '\r' << flush;  // clean last line of progress-output
        
        for(Host::Human &human : population->humans)
            human.flushReports();
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
