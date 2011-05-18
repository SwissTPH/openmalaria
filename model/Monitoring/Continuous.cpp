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
#include "util/errors.h"
#include "inputData.h"
#include "util/BoincWrapper.h"

#include <vector>
#include <map>
#include <fstream>
#include <boost/format.hpp>
#include <boost/math/nonfinite_num_facets.hpp>

namespace OM { namespace Monitoring {
    using namespace fastdelegate;
    using util::xml_scenario_error;
    
    /// Output file name
    string cts_filename;
    /// This is used to output some statistics in a tab-deliminated-value file.
    /// (It used to be csv, but German Excel can't open csv directly.)
    fstream ctsOStream;
    // Can't simply replace with ogzstream for compression: that doesn't support appending.
    
    /* Record last position in file (as position minus start), for checkpointing.
     * Don't use a streampos directly, because I'm not convinced we can save and
     * reload a streampos and use on a new file. */
    streamoff streamOff;
    streampos streamStart;
    
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
	ctsPeriod = ctsOpt.get().getPeriod();
	if( ctsPeriod < 1 )
	    throw xml_scenario_error("monitoring.continuous.period: must be >= 1 timestep");
	
        if( ctsOpt.get().getDuringInit().present() )
            duringInit = ctsOpt.get().getDuringInit().get();
        
        cts_filename = util::BoincWrapper::resolveFile("ctsout.txt");
        
	// This locale ensures uniform formatting of nans and infs on all platforms.
	locale old_locale;
	locale nfn_put_locale(old_locale, new boost::math::nonfinite_num_put<char>);
	ctsOStream.imbue( nfn_put_locale );
	ctsOStream.width (0);
	
	if( isCheckpoint ){
	    scnXml::OptionSet::OptionSequence sOSeq = ctsOpt.get().getOption();
	    for (scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
		map<string,Callback>::const_iterator reg_it = registered.find( it->getName() );
		if( reg_it == registered.end() )
		    throw xml_scenario_error( (boost::format("monitoring.continuous: no output \"%1%\"") %it->getName() ).str() );
		if( it->getValue() ){
		    toReport.push_back( reg_it->second.cb );
		}
	    }
	    
	    // When loading a check-point, we resume reporting to this file.
	    // Use "ate" mode and seek to desired pos.
	    ctsOStream.open (cts_filename.c_str(), ios::binary|ios::ate|ios::in|ios::out );
	    if( ctsOStream.fail() )
		throw util::checkpoint_error ("Continuous: resume error (no file)");
	    ctsOStream.seekp( 0, ios_base::beg );
	    streamStart = ctsOStream.tellp();
	    // we set position later, in staticCheckpoint
	}else{
	    if (util::BoincWrapper::fileExists(cts_filename.c_str())){
		// It could be from an old run. But we won't remove/truncate
		// existing files as a security precaution for running on BOINC.
		throw runtime_error ("File ctsout.txt exists!");
            }
	    
	    ctsOStream.open( cts_filename.c_str(), ios::binary|ios::out );
	    streamStart = ctsOStream.tellp();
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
	    ctsOStream << lineEnd << flush;
	    streamOff = ctsOStream.tellp() - streamStart;
	}
    }
    void Continuous::staticCheckpoint (ostream& stream){
        if( ctsPeriod == 0 )
            return;	// output disabled
	
	streamOff & stream;
    }
    void Continuous::staticCheckpoint (istream& stream){
        if( ctsPeriod == 0 )
            return;	// output disabled
	
	/* We attempt to resume output correctly after a reload by:
	 *
	 * (a) recording the last position, and relocating there
	 * (b) trying to avoid kills during writing of a line, by using
	 *  BOINC critical sections.
	 * 
	 * (Keeping results in memory until end of sim would be another,
	 * slightly safer, option, but loses real-time output.) */
	streamOff & stream;
	// We skip back to the last write-point, so anything written after the
	// last checkpoint will be repeated:
	ctsOStream.seekp( streamOff, ios_base::beg );
	
	if( ctsOStream.fail() )
	    throw util::checkpoint_error ("Continuous: resume error (bad pos/file)");
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
            return;	// output disabled
        if( !duringInit ){
            if( TimeStep::interventionPeriod < TimeStep(0) || TimeStep::interventionPeriod % ctsPeriod != 0 )
                return;
        } else {
            if( TimeStep::simulation % ctsPeriod != 0 )
                return;
            ctsOStream << TimeStep::simulation << '\t';
        }
	
	util::BoincWrapper::beginCriticalSection();	// see comment in staticCheckpoint
	
	ctsOStream << TimeStep::interventionPeriod;
	for( size_t i = 0; i < toReport.size(); ++i )
	    (toReport[i])( ctsOStream );
	// We must flush often to avoid temporarily outputting partial lines
	// (resulting in incorrect real-time graphs).
	ctsOStream << lineEnd << flush;
	
	streamOff = ctsOStream.tellp() - streamStart;
	util::BoincWrapper::endCriticalSection();
    }
} }
