/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2013 Swiss Tropical and Public Health Institute 
 * Copyright (C) 2005-2013 Liverpool School Of Tropical Medicine
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

#include "util/BoincWrapper.h"
#include "util/DocumentLoader.h"

#include "Global.h"
#include "Simulator.h"
#include "util/CommandLine.h"
#include "util/errors.h"

#include <cstdio>
#include <cerrno>

using namespace OM;

/** main() — initializes and shuts down BOINC, loads scenario XML and
 * runs simulation. */
int main(int argc, char* argv[]) {
    int exitStatus = EXIT_SUCCESS;
    string scenarioFile;
    
    try {
        util::set_gsl_handler();        // init
        
        scenarioFile = util::CommandLine::parse (argc, argv);   // parse arguments
        
        util::BoincWrapper::init();     // BOINC init
        
        // Load the scenario document:
        scenarioFile = util::CommandLine::lookupResource (scenarioFile);
        util::DocumentLoader documentLoader;
        util::Checksum cksum = documentLoader.loadDocument(scenarioFile);
        
        // Set up the simulator
        Simulator simulator( cksum, documentLoader.document() );
        
        // Save changes to the document if any occurred.
        documentLoader.saveDocument();
        
        // Write scenario checksum, only if simulation completed.
        cksum.writeToFile (util::BoincWrapper::resolveFile ("scenario.sum"));
        
        if ( !util::CommandLine::option(util::CommandLine::SKIP_SIMULATION) )
            simulator.start(documentLoader.document().getMonitoring());
        
        // We call boinc_finish before cleanup since it should help ensure
        // app isn't killed between writing output.txt and calling boinc_finish,
        // and may speed up exit.
        util::BoincWrapper::finish(exitStatus);	// Never returns
        
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
#ifdef WITHOUT_BOINC
        // Don't print this on BOINC, because if it's a problem we should find
        // it anyway!
        cerr << "This is likely an error in the C++ code. Please report!" << endl;
#endif
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
    
    // In case an exception was thrown, we call boinc_finish here:
    util::BoincWrapper::finish(exitStatus);	// Never returns
    return exitStatus;  // isn't actually reached, but avoids compiler warning
}
