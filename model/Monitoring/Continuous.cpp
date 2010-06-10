/* This file is part of OpenMalaria.
*
* Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Monitoring/Continuous.h"
#include "Monitoring/Survey.h"	// lineEnd
#include "util/errors.hpp"
#include "inputData.h"

#include <vector>
#include <map>
#include <fstream>
#include <boost/format.hpp>
#include <boost/math/nonfinite_num_facets.hpp>

namespace OM { namespace Monitoring {
    using namespace fastdelegate;
    using util::xml_scenario_error;
    
    /// Output file name
    const char* CTS_FILENAME = "ctsout.txt";
    /// This is used to output some statistics in a tab-deliminated-value file.
    /// (It used to be csv, but German Excel can't open csv directly.)
    ofstream ctsOStream;
    // Can't simply replace with ogzstream for compression: that doesn't support appending.
    
    // List of all registered callbacks (not used after init() runs)
    struct Callback {
	string titles;
	FastDelegate1<ostream&> cb;
    };
    map<string,Callback> registered;
    
    // List that we report.
    vector< FastDelegate1<ostream&> > toReport;
    int ctsPeriod = 0;
    bool duringInit = false;
    
    /* Initialise: enable outputs registered and requested in XML.
     * Search for Continuous::registerCallback to see outputs available. */
    void Continuous::init (bool isCheckpoint) {
	const scnXml::Monitoring::ContinuousOptional& ctsOpt =
	    InputData().getMonitoring().getContinuous();
	if( ctsOpt.present() == false ) {
	    ctsPeriod = 0;
	    return;
	}
	ctsPeriod = ctsOpt.get().getPeriod() / Global::interval;
	if( ctsPeriod < 1 )
	    throw xml_scenario_error("monitoring.continuous.period: must be > 1 timestep");
	else if( ctsOpt.get().getPeriod() % Global::interval != 0 )
	    cerr << "Warning: monitoring.continuous.period not a whole number of timesteps" << endl;
	
        if( ctsOpt.get().getDuringInit().present() )
            duringInit = ctsOpt.get().getDuringInit().get();
        
	// This locale ensures uniform formatting of nans and infs on all platforms.
	locale old_locale;
	locale nfn_put_locale(old_locale, new boost::math::nonfinite_num_put<char>);
	ctsOStream.imbue( nfn_put_locale );
	ctsOStream.width (0);
	
	if( isCheckpoint )
	    // When loading a check-point, we resume reporting to this file.
	    // Don't worry about writing to an existing file (like for "output.txt"): we're appending.
	    ctsOStream.open (CTS_FILENAME, ios::app | ios::binary);
	else{
	    ifstream test (CTS_FILENAME);
	    if (test.is_open())
		// It could be from an old run. But we won't remove/truncate
		// existing files as a security precaution for running on BOINC.
		throw runtime_error ("File vector.txt exists!");
	    
	    ctsOStream.open( CTS_FILENAME );
	    ctsOStream << "##\t##" << endl;	// live-graph needs a deliminator specifier when it's not a comma
	    
	    if( duringInit )
                ctsOStream << "simulation time\t";
	    ctsOStream << "timestep";
	    scnXml::OptionSet::OptionSequence sOSeq = ctsOpt.get().getOption();
	    for (scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
		map<string,Callback>::const_iterator reg_it = registered.find( it->getName() );
		if( reg_it == registered.end() )
		    throw xml_scenario_error( (boost::format("monitoring.continuous: no output \"%1%\"") %it->getName() ).str() );
		if( it->getValue() ){
		    ctsOStream << reg_it->second.titles;
		    toReport.push_back( reg_it->second.cb );
		}
	    }
	    ctsOStream << lineEnd;
	}
    }
    
	void Continuous::registerCallback (string optName, string titles, fastdelegate::FastDelegate1<ostream&> outputCb){
	Callback s;
	s.titles = titles;
	s.cb = outputCb;
	assert(registered.count(optName) == 0);	// name clash/registered twice?
	registered[optName] = s;
    }
    
    void Continuous::update (){
        if( ctsPeriod == 0 )
            return;
        if( !duringInit ){
            if( Global::timeStep < 0 || Global::timeStep % ctsPeriod != 0 )
                return;
        } else {
            if( Global::simulationTime % ctsPeriod != 0 )
                return;
            ctsOStream << Global::simulationTime << '\t';
        }
	
	ctsOStream << Global::timeStep;
	for( size_t i = 0; i < toReport.size(); ++i )
	    (toReport[i])( ctsOStream );
	ctsOStream << lineEnd;
    }
} }
