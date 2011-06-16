/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include "inputData.h"	//Include parser for the input

#include "Global.h"
#include "Simulation.h"
#include "util/CommandLine.h"

#include <cstdio>
#include <cerrno>

/** main() — initializes and shuts down BOINC, loads scenario XML and
 * runs simulation. */
int main(int argc, char* argv[]) {
    int exitStatus = EXIT_SUCCESS;
    
    try {
        string scenario_name =
            OM::util::CommandLine::parse (argc, argv);
	
        OM::util::BoincWrapper::init();
	
        scenario_name = OM::util::CommandLine::lookupResource (scenario_name);
        OM::util::Checksum cksum = OM::InputData.createDocument(scenario_name);
	
	OM::Simulation simulation (cksum);	// constructor runs; various initialisations
	
	// Save changes to the document if any occurred.
        OM::InputData.saveDocument();
	
	if ( !OM::util::CommandLine::option(OM::util::CommandLine::SKIP_SIMULATION) )
	    simulation.start();
	
	// We call boinc_finish before cleanup since it should help ensure
	// app isn't killed between writing output.txt and calling boinc_finish,
	// and may speed up exit.
	OM::util::BoincWrapper::finish(exitStatus);	// Never returns
	
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
        cerr << "XSD Exception: " << e.what() << '\n' << e << endl;
        exitStatus = OM::util::Error::XSD;
    } catch (const OM::util::checkpoint_error& e) {
        cerr << "Checkpoint exception: " << e.what() << endl;
        cerr << e << flush;
        exitStatus = e.getCode();
    } catch (const OM::util::traced_exception& e) {
        cerr << "Exception: " << e.what() << endl;
        cerr << e << flush;
        exitStatus = e.getCode();
    } catch (const OM::util::xml_scenario_error& e) {
        cerr << "Error in scenario XML file: " << e.what() << endl;
        exitStatus = e.getCode();
    } catch (const OM::util::base_exception& e) {
        cerr << "Exception: " << e.what() << endl;
        exitStatus = e.getCode();
    } catch (const exception& e) {
        cerr << "Exception: " << e.what() << endl;
        exitStatus = EXIT_FAILURE;
    } catch (...) {
        cerr << "Unknown exception" << endl;
        exitStatus = EXIT_FAILURE;
    }
    
    // If we get to here, we already know an error occurred.
    if( errno != 0 )
	std::perror( "OpenMalaria" );
    
    // Free memory (though usually we don't bother at exit to save time)
    OM::InputData.freeDocument();
    
    // In case an exception was thrown, we call boinc_finish here:
    OM::util::BoincWrapper::finish(exitStatus);	// Never returns
}
