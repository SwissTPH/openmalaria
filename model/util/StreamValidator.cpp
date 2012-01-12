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
#include "util/CommandLine.h"
#include "util/errors.h"

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
    // Versions of checkpoint operator for deque, but with increased max
    // length. They also write length as a 32-bit unsigned integer; in theory
    // size_t could be used for checkpointing but that's not really necessary.
    template<class T>
    void operator& (deque<T> x, ostream& stream) {
        // Write size as a uint32_t for cross-platform compatibility; first
        // check size isn't greater than what a 32-bit uint can represent.
        uint32_t size32 = x.size();
        size_t size_check = size32;
        if( size_check != x.size() )
            throw traced_exception( "stream too long!" );
        size32 & stream;       // write 32-bit unsigned int (size_t is platform dependent!)
        BOOST_FOREACH (T& y, x) {
            y & stream;
        }
    }
    template<class T>
    void operator& (deque<T>& x, istream& stream) {
        uint32_t size32;
        size32 & stream; // read uint32_t, not whatever size_t is
        validateListSize (size32, 100000000L);
        x.resize( size32 );
        BOOST_FOREACH (T& y, x) {
            y & stream;
        }
    }
}

void StreamValidatorType::saveStream() {
    if( storeMode ){
	ofstream f_str( OM_SV_FILE, ios::out | ios::binary );
	if( !f_str.is_open() )
	    throw traced_exception( "unable to write " OM_SV_FILE, Error::FileIO );
	f_str.write( reinterpret_cast<const char*>(&OM_SV_HEAD), sizeof(char)*4 );
        
        stream & f_str;
        
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
	throw traced_exception( (boost::format("unable to read %1%") %file).str(), Error::FileIO );
    char head[4];
    f_str.read( reinterpret_cast<char*>(&head), sizeof(char)*4 );
    if( memcmp( &OM_SV_HEAD, &head, sizeof(char)*4 ) != 0 )
	throw traced_exception( (boost::format("%1% is not a valid StreamValidator file") %file).str(), Error::FileIO );
    
    stream & f_str;
    
    f_str.ignore (numeric_limits<streamsize>::max()-1);	// skip to end of file
    if (f_str.gcount () != 0) {
	ostringstream msg;
	msg << file<<" has " << f_str.gcount() << " bytes remaining." << endl;
	throw traced_exception (msg.str(), Error::FileIO);
    } else if (f_str.fail())
	throw traced_exception ("StreamValidator load error", Error::FileIO);
    
    f_str.close();
    readIt = stream.begin();
}

void StreamValidatorType::handle( SVType value ){
    if( storeMode ){
	stream.push_back( value );
    }else{
	if( value != *readIt ){
	    // Attach a debugger with a breakpoint here to get the backtrace,
	    // if the one caught by traced_exception isn't enough.
	    throw traced_exception ("StreamValidator: out of sync!");
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

// ———  Our cross-platform consistent-result hasing functions  ——
namespace CPCH {
    SVType toSVType(boost::uint32_t x){
        return static_cast<SVType>(x);
    }
    SVType toSVType(boost::int32_t x){
        return static_cast<SVType>(x);
    }
    SVType toSVType(boost::uint64_t x){
        return static_cast<SVType>(x);
    }
    SVType toSVType(boost::int64_t x){
        return static_cast<SVType>(x);
    }
    //NOTE: hashing of floats and doubles isn't very good or fast:
    SVType toSVType(float x){
        return x*1e4f;
    }
    SVType toSVType(double x){
        return x*1e4;
    }
}


} }
#endif
