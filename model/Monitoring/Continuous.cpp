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

// must be included first to avoid compiler error (incompatible with fpclassify which is included from Global.h?)
#include <boost/math/nonfinite_num_facets.hpp>

#include "Monitoring/Continuous.h"
#include "mon/info.h"   // lineEnd
#include "util/errors.h"
#include "util/BoincWrapper.h"
#include "util/CommandLine.h"
#include "schema/monitoring.h"

#include <vector>
#include <map>
#include <fstream>
#include <boost/format.hpp>
#include <gzstream/gzstream.h>

namespace OM { namespace Monitoring {
    using namespace fastdelegate;
    using util::xml_scenario_error;
    
    /// File we send uncompressed output to
    string cts_filename;
#ifndef WITHOUT_BOINC
    /// At end of simulation, compress the output file. Don't compress data as
    /// it's output, because gzstream doesn't support seeking, which is needed
    /// for checkpoint resume.
    string compressedCtsoutName;
#endif
    /// This is used to output some statistics in a tab-deliminated-value file.
    /// (It used to be csv, but German Excel can't open csv directly.)
    fstream ctsOStream;
    
    /* Record last position in file (as position minus start), for checkpointing.
     * Don't use a streampos directly, because I'm not convinced we can save and
     * reload a streampos and use on a new file. */
    streamoff streamOff;
    streampos streamStart;
    
    // List of all registered callbacks (not used after init() runs)
    class Callback {
    protected:
        Callback( const string& t ) : titles(t) {}
    public:
        virtual ~Callback() {}
        string titles;
        virtual void call( const Population&, ostream& ) =0;
    };
    class Callback1 : public Callback {
	FastDelegate1<ostream&> cb;
    public:
        Callback1( const string& t, FastDelegate1<ostream&> outputCb ) :
            Callback(t), cb( outputCb ) {}
        virtual void call( const Population&, ostream& stream ){
            cb( stream );
        }
    };
    class Callback2Pop : public Callback {
        FastDelegate2<const Population&,ostream&> cb;
    public:
        Callback2Pop( const string& t, FastDelegate2<const Population&,ostream&> outputCb ) :
            Callback(t), cb( outputCb ) {}
        virtual void call( const Population& pop, ostream& stream ){
            cb( pop, stream );
        }
    };
    typedef map<string,Callback*> registered_t;
    registered_t registered;
    
    // List that we report.
    vector< Callback* > toReport;
    SimTime ctsPeriod = sim::zero();
    bool duringInit = false;
    
    
    ContinuousType Continuous;
    
    ContinuousType::~ContinuousType (){
        // free memory
        toReport.clear();
        for( registered_t::iterator it = registered.begin(); it != registered.end(); ++it )
            delete it->second;
   }
   
    /* Initialise: enable outputs registered and requested in XML.
     * Search for Continuous::registerCallback to see outputs available. */
    void ContinuousType::init (const scnXml::Monitoring& monitoring, bool isCheckpoint) {
	const scnXml::Monitoring::ContinuousOptional& ctsOpt = monitoring.getContinuous();
	if( ctsOpt.present() == false ) {
	    ctsPeriod = sim::zero();
	    return;
	}
	try{
            //NOTE: if changing XSD, this should not have a default unit:
            ctsPeriod = UnitParse::readShortDuration( ctsOpt.get().getPeriod(), UnitParse::STEPS );
            if( ctsPeriod < sim::oneTS() )
                throw util::format_error("must be >= 1 time step");
        }catch( const util::format_error& e ){
            throw xml_scenario_error( string("monitoring/continuous/period: ").append(e.message()) );
        }
	
        if( ctsOpt.get().getDuringInit().present() )
            duringInit = ctsOpt.get().getDuringInit().get();
        
        cts_filename = util::BoincWrapper::resolveFile(util::CommandLine::getCtsoutName());
#ifndef WITHOUT_BOINC
        // redirect output to a temporary file; copy and compress this to the
        // final output at end of simulation
        compressedCtsoutName = cts_filename;
        cts_filename = "ctsout_temp.txt";
#endif
        
	// This locale ensures uniform formatting of nans and infs on all platforms.
	locale old_locale;
	locale nfn_put_locale(old_locale, new boost::math::nonfinite_num_put<char>);
	ctsOStream.imbue( nfn_put_locale );
	ctsOStream.width (0);
	
	if( isCheckpoint ){
	    scnXml::OptionSet::OptionSequence sOSeq = ctsOpt.get().getOption();
	    for (scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
		registered_t::const_iterator reg_it = registered.find( it->getName() );
		if( reg_it == registered.end() )
		    throw xml_scenario_error( (boost::format("monitoring.continuous: no output \"%1%\"") %it->getName() ).str() );
		if( it->getValue() ){
		    toReport.push_back( reg_it->second );
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
#ifndef WITHOUT_BOINC
	    if (util::BoincWrapper::fileExists(cts_filename.c_str())){
		// It could be from an old run. But we won't remove/truncate
		// existing files as a security precaution for running on BOINC.
		throw util::base_exception (string("File ").append(cts_filename).append(" exists!"),util::Error::FileExists);
            }
#endif
	    
	    ctsOStream.open( cts_filename.c_str(), ios::binary|ios::out );
	    streamStart = ctsOStream.tellp();
	    ctsOStream << "##\t##" << endl;	// live-graph needs a deliminator specifier when it's not a comma
	    
	    if( duringInit )
                ctsOStream << "simulation time\t";
	    ctsOStream << "timestep";   //TODO: change to days or remove or leave?
	    scnXml::OptionSet::OptionSequence sOSeq = ctsOpt.get().getOption();
	    for (scnXml::OptionSet::OptionConstIterator it = sOSeq.begin(); it != sOSeq.end(); ++it) {
		registered_t::const_iterator reg_it = registered.find( it->getName() );
		if( reg_it == registered.end() )
		    throw xml_scenario_error( (boost::format("monitoring.continuous: no output \"%1%\"") %it->getName() ).str() );
		if( it->getValue() ){
		    ctsOStream << reg_it->second->titles;
		    toReport.push_back( reg_it->second );
		}
	    }
	    ctsOStream << mon::lineEnd << flush;
	    streamOff = ctsOStream.tellp() - streamStart;
	}
    }
   
   void ContinuousType::finalise() {
         if( ctsPeriod == sim::zero() )
             return;     // output disabled
#ifndef WITHOUT_BOINC
        if (util::BoincWrapper::fileExists(compressedCtsoutName.c_str())){
            throw util::base_exception(string("File ").append(compressedCtsoutName).append(" exists!"),util::Error::FileExists);
        }
        ctsOStream.close();
        ifstream origFile(cts_filename.c_str());
        if( !origFile.is_open() ){
            throw util::base_exception(string("Temporary file ").append(cts_filename).append(" not found!"),util::Error::FileIO);
        }
        ogzstream finalFile(compressedCtsoutName.c_str());
        finalFile << origFile.rdbuf();
#endif
    }
    void ContinuousType::checkpoint (ostream& stream){
        if( ctsPeriod == sim::zero() )
            return;	// output disabled
	
	streamOff & stream;
    }
    void ContinuousType::checkpoint (istream& stream){
        if( ctsPeriod == sim::zero() )
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
    
    void ContinuousType::registerCallback (string optName, string titles,
            fastdelegate::FastDelegate1<ostream&> outputCb){
        assert(registered.count(optName) == 0); // name clash/registered twice?
	registered[optName] = new Callback1( titles, outputCb );
    }
    void ContinuousType::registerCallback (string optName, string titles,
            fastdelegate::FastDelegate2<const Population&,ostream&> outputCb){
        assert(registered.count(optName) == 0); // name clash/registered twice?
        registered[optName] = new Callback2Pop( titles, outputCb );
    }
    
    void ContinuousType::update (const Population& population){
        if( ctsPeriod == sim::zero() )
            return;	// output disabled
        if( !duringInit ){
            if( sim::intervNow() < sim::zero() || mod_nn(sim::intervNow(), ctsPeriod) != sim::zero() )
                return;
        } else {
            if( mod_nn(sim::now(), ctsPeriod) != sim::zero() )
                return;
            ctsOStream << sim::now().inSteps() << '\t';
        }
	
	util::BoincWrapper::beginCriticalSection();	// see comment in staticCheckpoint
	
        if( duringInit && sim::intervNow() < sim::zero() ){
            ctsOStream << "nan";
        }else{
            ctsOStream << sim::intervNow().inSteps();
        }
	for( size_t i = 0; i < toReport.size(); ++i )
	    toReport[i]->call( population, ctsOStream );
	// We must flush often to avoid temporarily outputting partial lines
	// (resulting in incorrect real-time graphs).
	ctsOStream << mon::lineEnd << flush;
	
	streamOff = ctsOStream.tellp() - streamStart;
	util::BoincWrapper::endCriticalSection();
    }
} }
