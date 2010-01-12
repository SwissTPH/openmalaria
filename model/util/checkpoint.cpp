/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "util/checkpoint.hpp"
#include "util/errors.hpp"

#include <limits>
#include <sstream>
using namespace std;

namespace OM { namespace util { namespace checkpoint {
    
    // Header test constants
    const unsigned int h_BOM = 0x50434D4F;	// "OMCP" in little-endian: OpenMalaria CheckPoint
    const bool h_b = true;
    const unsigned char h_c = 0xA5;	// binary: 10100101; don't care if char is signed as long as first bit is read and written correctly
    const double h_n0 = -0.0;
    const double h_nan = numeric_limits<double>::quiet_NaN();
    
    
    void validateListSize (long length) {
	if (length < 0 || length > 1000) {
	    ostringstream s;
	    s << "List length out of range: " << length;
	    throw checkpoint_error(s.str());
	}
    }
    
    /** @brief Binary checkpointing
     *
     * Writes value in direct binary (not suitible for pointers). Currently
     * machine dependent (don't transfer a checkpoint to another machine).
     */
    //@{
    template<class T>
    void binary_write (T x, ostream& stream) {
	stream.write (reinterpret_cast<char*>(&x), sizeof(x));
    }
    template<class T>
    void binary_read (T& x, istream& stream) {
	stream.read (reinterpret_cast<char*>(&x), sizeof(x));
	if (!stream || stream.gcount() != sizeof(x))
	    throw checkpoint_error ("stream read error binary");
    }
    //@}
    
    void staticChecks () {
	BOOST_STATIC_ASSERT (sizeof(char) == 1);
	BOOST_STATIC_ASSERT (sizeof(short int) == 2);
	BOOST_STATIC_ASSERT (sizeof(int) == 4);
	BOOST_STATIC_ASSERT (sizeof(float) == 4);
	BOOST_STATIC_ASSERT (sizeof(double) == 8);
    }
    
    void header (ostream& stream) {
	staticChecks();
	
	binary_write (h_BOM, stream);
	binary_write (h_b, stream);
	binary_write (h_c, stream);
	binary_write (h_n0, stream);
	binary_write (h_nan, stream);
    }
    void header (istream& stream) {
	staticChecks ();
	unsigned int BOM;
	bool b;
	unsigned char c;
	double n0;
	double nan;
	
	binary_read (BOM, stream);
	binary_read (b, stream);
	binary_read (c, stream);
	binary_read (n0, stream);
	binary_read (nan, stream);
	
	// Check. Use binary check for doubles since it's not the same as numeric ==
	if (BOM != h_BOM ||
	    b != h_b ||
	    c != h_c ||
	    memcmp (&n0, &h_n0, sizeof(double)) ||
	    memcmp (&nan, &h_nan, sizeof(double))
	    )
	    throw checkpoint_error ("invalid header");
    }
    
    ///@brief Operator& for simple data-types
    //@{
    void operator& (bool x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (bool& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (signed char x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (signed char& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (short x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (short& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (int x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (int& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (long x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (long& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (long long x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (long long& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (unsigned char x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (unsigned char& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (unsigned short x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (unsigned short& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (unsigned int x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (unsigned int& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (unsigned long x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (unsigned long& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (unsigned long long x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (unsigned long long& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (float x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (float& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (double x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (double& x, istream& stream) {
	binary_read (x, stream);
    }
    
    void operator& (long double x, ostream& stream) {
	binary_write (x, stream);
    }
    void operator& (long double& x, istream& stream) {
	binary_read (x, stream);
    }
    //@}
    
    // string
    void operator& (string x, ostream& stream) {
	x.length() & stream;
	stream.write (x.c_str(), x.length());
    }
    void operator& (string& x, istream& stream) {
	size_t len;
	len & stream;
	x.resize (len);
	stream.read (&x[0], x.length());
	if (!stream || stream.gcount() != streamsize(len))
	    throw checkpoint_error ("stream read error string");
    }
    
    // map<S,T> — template doesn't work on gcc
    void operator& (map<string,double> x, ostream& stream) {
	x.size() & stream;
	for (map<string,double>::const_iterator pos = x.begin (); pos != x.end() ; ++pos) {
	    pos->first & stream;
	    pos->second & stream;
	}
    }
    void operator& (map<string,double>& x, istream& stream) {
	size_t l;
	l & stream;
	validateListSize (l);
	x.clear ();
	map<string,double>::iterator pos = x.begin ();
	for (size_t i = 0; i < l; ++i) {
	    string s;
	    double t;
	    s & stream;
	    t & stream;
	    pos = x.insert (pos, make_pair (s,t));
	}
    }
    
} } }
