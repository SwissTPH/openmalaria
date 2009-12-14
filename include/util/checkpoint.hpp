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
#ifndef OM_util_checkpoint
#define OM_util_checkpoint

#include "util/errors.hpp"
#include <iostream>
#include <boost/foreach.hpp>
#include <list>
using namespace std;

/** @brief Checkpointing utility functions
 *
 * These functions are intended to facilitate checkpointing.
 * 
 * My current rule-of-thumb is to checkpoint all non-static data but not static data in general.
 * 
 * Data should be checkpointed with the & operator. Non-derived classes should implement a
 * checkpointing function like the following:
 * 
 * \code
/// Checkpointing
template<class S>
void operator& (S& stream) {
    var1 & stream;
    var2 & stream;
}
 * \endcode
 * 
 * (This template will be implemented for istream and ostream types.)
 * Derived classes should implement the same templated function, but calling a virtual checkpoint
 * function defined as the following:
 * 
 * \code
virtual checkpoint (istream& stream);
virtual checkpoint (ostream& stream);
 * \endcode
 *
 * Although this needs to be implemented twice, in most cases one implementation will just be a copy
 * of the other with istream replaced with ostream as the argument type.
 * 
 * One should declare "using namespace OM::util::checkpoint" to use the operator& functions defined
 * here. */
namespace OM { namespace util { namespace checkpoint {
    
    ///@brief Utility functions
    //@{
  /** For checking the number of elements in a list is sensible when loading
   * checkpoints.
   *
   * While loading checkpoints, if memory is allocated in a loop using a
   * variable loaded from the checkpoint file, it should be checked that length
   * is sensible. This function is a convenient way to do that (although the
   * limit isn't large enough to account for population size).
   * 
   * The reason for this check is that if the checkpoint is read wrongly,
   * variables often get completely wrong values which can cause memory
   * allocation to grind the computer to a halt. In this case the values are
   * usually bad enough that a lenient check will still catch them. */
    void validateListSize (long length);
    //@}
    
    ///@brief Operator& for simple data-types
    //@{
    inline void operator& (bool x, ostream& stream) { stream << x << endl; }
    inline void operator& (bool& x, istream& stream) { stream >> x; }
    
    inline void operator& (signed char x, ostream& stream) { stream << static_cast<short>(x) << endl; }
    inline void operator& (signed char& x, istream& stream) {
	short s;
	stream >> s;
	x = s;
    }
    
    inline void operator& (short x, ostream& stream) { stream << x << endl; }
    inline void operator& (short& x, istream& stream) { stream >> x; }
    
    inline void operator& (int x, ostream& stream) { stream << x << endl; }
    inline void operator& (int& x, istream& stream) { stream >> x; }
    
    inline void operator& (long x, ostream& stream) { stream << x << endl; }
    inline void operator& (long& x, istream& stream) { stream >> x; }
    
    inline void operator& (long long x, ostream& stream) { stream << x << endl; }
    inline void operator& (long long& x, istream& stream) { stream >> x; }
    
    inline void operator& (unsigned char x, ostream& stream) { stream << static_cast<unsigned short>(x) << endl; }
    inline void operator& (unsigned char& x, istream& stream) {
	unsigned short us;
	stream >> us;
	x = us;
    }
    
    inline void operator& (unsigned short x, ostream& stream) { stream << x << endl; }
    inline void operator& (unsigned short& x, istream& stream) { stream >> x; }
    
    inline void operator& (unsigned int x, ostream& stream) { stream << x << endl; }
    inline void operator& (unsigned int& x, istream& stream) { stream >> x; }
    
    inline void operator& (unsigned long x, ostream& stream) { stream << x << endl; }
    inline void operator& (unsigned long& x, istream& stream) { stream >> x; }
    
    inline void operator& (unsigned long long x, ostream& stream) { stream << x << endl; }
    inline void operator& (unsigned long long& x, istream& stream) { stream >> x; }
    
    inline void operator& (float x, ostream& stream) { stream << x << endl; }
    inline void operator& (float& x, istream& stream) { stream >> x; }
    
    inline void operator& (double x, ostream& stream) { stream << x << endl; }
    inline void operator& (double& x, istream& stream) { stream >> x; }
    
    inline void operator& (long double x, ostream& stream) { stream << x << endl; }
    inline void operator& (long double& x, istream& stream) { stream >> x; }
    //@}
    
    /** @brief Operator& for pointers
     *
     * These could be implemented for pointers but require reading objects in the opposite order to
     * writing them.
     * 
     * Matching up pointers with addresses of objects not stored in references could be more
     * challanging though.
     * 
     * Can't really deal with referenced objects properly. */
    //@{
    //@}
    
    /** @brief Operator& for std::string */
    //@{
    void operator& (string x, ostream& stream);
    void operator& (string& x, istream& stream);
    //@}
    
    ///@brief Operator& for stl containers
    //@{
    template<class T>
    void operator& (vector<T> x, ostream& stream) {
	x.size() & stream;
	BOOST_FOREACH (T& y, x) {
	    y & stream;
	}
    }
    template<class T>
    void operator& (vector<T>& x, istream& stream) {
	size_t l;
	l & stream;
	validateListSize (l);
	x.resize (l);
	BOOST_FOREACH (T& y, x) {
	    y & stream;
	}
    }
    
    template<class T>
    void operator& (list<T> x, ostream& stream) {
	x.size() & stream;
	BOOST_FOREACH (T& y, x) {
	    y & stream;
	}
    }
    template<class T>
    void operator& (list<T>& x, istream& stream) {
	size_t l;
	l & stream;
	validateListSize (l);
	x.resize (l);
	BOOST_FOREACH (T& y, x) {
	    y & stream;
	}
    }
    //@}
    
} } }	// end of namespaces
#endif
