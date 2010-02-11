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

#include <iostream>
#include <boost/foreach.hpp>
#include <vector>
#include <list>
#include <map>
using namespace std;

/** @brief Checkpointing utility functions
 *
 * These functions are intended to facilitate checkpointing. They were inspired
 * by the boost::serialization library (which proved to be easier not to use,
 * due to the additional library dependency).
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
    /** Perform important checks on checkpoint file format.
     *
     * Call at beginning of read/write. */
    void header (ostream& stream);
    void header (istream& stream);	///< ditto
    
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
    void operator& (bool x, ostream& stream);
    void operator& (bool& x, istream& stream);
    
    void operator& (signed char x, ostream& stream);
    void operator& (signed char& x, istream& stream);
    
    void operator& (short x, ostream& stream);
    void operator& (short& x, istream& stream);
    
    void operator& (int x, ostream& stream);
    void operator& (int& x, istream& stream);
    
    void operator& (long x, ostream& stream);
    void operator& (long& x, istream& stream);
    
    void operator& (long long x, ostream& stream);
    void operator& (long long& x, istream& stream);
    
    void operator& (unsigned char x, ostream& stream);
    void operator& (unsigned char& x, istream& stream);
    
    void operator& (unsigned short x, ostream& stream);
    void operator& (unsigned short& x, istream& stream);
    
    void operator& (unsigned int x, ostream& stream);
    void operator& (unsigned int& x, istream& stream);
    
    void operator& (unsigned long x, ostream& stream);
    void operator& (unsigned long& x, istream& stream);
    
    void operator& (unsigned long long x, ostream& stream);
    void operator& (unsigned long long& x, istream& stream);
    
    void operator& (float x, ostream& stream);
    void operator& (float& x, istream& stream);
    
    void operator& (double x, ostream& stream);
    void operator& (double& x, istream& stream);
    
    void operator& (long double x, ostream& stream);
    void operator& (long double& x, istream& stream);
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
    /// Version of above taking an element to initialize each element from.
    template<class T>
    void checkpoint (vector<T>& x, istream& stream, T templateInstance) {
	size_t l;
	l & stream;
	validateListSize (l);
	x.resize (l, templateInstance);
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
    
    /* The following templates don't work on gcc, so you'll need to provide them as non-templates.
    template<class S, class T>
    void operator& (map<S,T> x, ostream& stream) {
	x.size() & stream;
	BOOST_FOREACH (T& y, x) {
	    y->first & stream;
	    y->second & stream;
	}
    }
    template<class S, class T>
    void operator& (map<S,T> x, istream& stream) {
	size_t l;
	l & stream;
	validateListSize (l);
	x.clear ();
	map<S,T>::iterator pos = x.begin ();
	for (size_t i = 0; i < l; ++i) {
	    S s;
	    T t;
	    s & stream;
	    t & stream;
	    pos = x.insert (pos, make_pair (s,t));
	}
    }
    */
    
    void operator& (map<string,double> x, ostream& stream);
    void operator& (map<string, double >& x, istream& stream);
    
    void operator& (multimap<double,double> x, ostream& stream);
    void operator& (multimap<double, double>& x, istream& stream);
    //@}
    
} } }	// end of namespaces
#endif
