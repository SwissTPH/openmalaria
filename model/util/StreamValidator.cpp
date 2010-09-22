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

#include "util/StreamValidator.h"
#include "util/CommandLine.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <boost/format.hpp>

// Compile-time optional
#ifdef OM_STREAM_VALIDATOR
namespace OM { namespace util {

#define OM_SV_FILE "StreamValidator"
const char OM_SV_HEAD[5] = "OMSV";

// These operators overload against others in the checkpoint namespace.
// Putting these in the same namespace is an easy solution, if not very standard.
namespace checkpoint {
    // versions of checkpoint operator for deque,
    // but with increased max length
    template<class T>
    void operator& (deque<T> x, ostream& stream) {
	x.size() & stream;
	BOOST_FOREACH (T& y, x) {
	    y & stream;
	}
    }
    template<class T>
    void operator& (deque<T>& x, istream& stream) {
	size_t l;
	l & stream;
	validateListSize (l, 1e8);
	x.resize( l );
	BOOST_FOREACH (T& y, x) {
	    y & stream;
	}
    }
}

void StreamValidatorType::saveStream() {
    if( storeMode ){
	ofstream f_str( OM_SV_FILE, ios::out | ios::binary );
	if( !f_str.is_open() )
	    throw runtime_error( "unable to write " OM_SV_FILE );
	f_str.write( reinterpret_cast<const char*>(&OM_SV_HEAD), sizeof(char)*4 );
	stream & f_str;	// use checkpointing operation
	f_str.close();
    }else{
	if( readIt != stream.end() )
	    cerr << "StreamValidator: not at end!" << endl;
    }
}

void StreamValidatorType::loadStream( const string& path ){
    string file = CommandLine::lookupResource( path );
    storeMode = false;
    ifstream f_str( file.c_str(), ios::in | ios::binary );
    if( !f_str.is_open() )
	throw runtime_error( (boost::format("unable to read %1%") %file).str() );
    char head[4];
    f_str.read( reinterpret_cast<char*>(&head), sizeof(char)*4 );
    if( memcmp( &OM_SV_HEAD, &head, sizeof(char)*4 ) != 0 )
	throw runtime_error( (boost::format("%1% is not a valid StreamValidator file") %file).str() );
    stream & f_str;	// checkpointing operator
    
    f_str.ignore (numeric_limits<streamsize>::max()-1);	// skip to end of file
    if (f_str.gcount () != 0) {
	ostringstream msg;
	msg << file<<" has " << f_str.gcount() << " bytes remaining." << endl;
	throw runtime_error (msg.str());
    } else if (f_str.fail())
	throw runtime_error ("StreamValidator load error");
    
    f_str.close();
    readIt = stream.begin();
}

void StreamValidatorType::handle( size_t value ){
    if( storeMode ){
	stream.push_back( value );
    }else{
	if( value != *readIt ){
	    // Attach a debugger with a breakpoint here to get the backtrace.
	    cerr << "StreamValidator: out of sync!" << endl;
	}
	++readIt;
    }
}

void StreamValidatorType::checkpoint ( istream& cp_str ){
    storeMode & cp_str;
    stream & cp_str;
    size_t dist;
    dist & cp_str;
    readIt = stream.begin() + dist;
}
void StreamValidatorType::checkpoint ( ostream& cp_str ) const{
    storeMode & cp_str;
    stream & cp_str;
    size_t dist = readIt - stream.begin();
    dist & cp_str;
}

StreamValidatorType StreamValidator;

} }
#endif
