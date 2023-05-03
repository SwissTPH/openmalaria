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

#include "util/errors.h"
#include "util/CommandLine.h"
#include "util/ModelOptions.h"
#include "util/StreamValidator.h"
#include "util/DocumentLoader.h"

#include "mon/Continuous.h"
#include "mon/management.h"

#include "interventions/InterventionManager.hpp"
#include "Clinical/ClinicalModel.h"

#include "checkpoint.h"

#include "schema/scenario.h"

#include "Model.h"

#include <cerrno>

namespace OM {
    using mon::Continuous;
    using interventions::InterventionManager;
    using Transmission::TransmissionModel;
}

using namespace OM;

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

        std::unique_ptr<Model> model = unique_ptr<Model>(model::create(*scenario));
        
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
        
        estEndTime = model->humanWarmupLength + model->transmission->expectedInitDuration() + (sim::endDate() - sim::startDate()) + sim::oneTS();
        assert( estEndTime + sim::never() < sim::zero() );

        if (startedFromCheckpoint)
        {
            Continuous.init(scenario->getMonitoring(), true);
            checkpoint::read(checkpointFileName, endTime, estEndTime, *model->population, *model->transmission);
        }
        else
        {
            Continuous.init(scenario->getMonitoring(), false);
            model->population->createInitialHumans();
            model->transmission->init2(model->population->humans);
        
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
            endTime = model->humanWarmupLength;
            model::run(*model, endTime, estEndTime, "Warmup");

            /** Transmission init phase:
             * Fit the emergence rate to the input EIR */
            SimTime iterate = model->transmission->initIterate();
            while(iterate > sim::zero())
            {
                endTime = endTime + iterate;
                // adjust estimation of final time step: end of current period + length of main phase
                estEndTime = endTime + (sim::endDate() - sim::startDate()) + sim::oneTS();
                model::run(*model, endTime, estEndTime, "EIR Calibration");
                iterate = model->transmission->initIterate();
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
            model->population->recentBirths = 0;
            model->transmission->summarize(); // Only to reset TransmissionModel::inoculationsPerAgeGroup
            mon::initMainSim();

            if(util::CommandLine::option (util::CommandLine::CHECKPOINT))
            {
                checkpoint::write(startedFromCheckpoint, checkpointFileName, endTime, estEndTime, *model->population, *model->transmission);
                if( util::CommandLine::option (util::CommandLine::CHECKPOINT_STOP) )
                    throw util::cmd_exception ("Checkpoint test: checkpoint written", util::Error::None);
            }
        }

        // Main phase loop
        model::run(*model, endTime, estEndTime, "Intervention period");
       
        cerr << '\r' << flush;  // clean last line of progress-output
        
        for(Host::Human &human : model->population->humans)
            human.clinicalModel->flushReports();
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
