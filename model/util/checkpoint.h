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

#ifndef OM_util_checkpoint
#define OM_util_checkpoint

#ifndef Hmod_Global
#error "Please include Global.h not checkpoint.h directly."
// otherwise "using ..." declaration in Global.h won't work
#endif

#include <iostream>
#include <boost/foreach.hpp>
#include <vector>
#include <list>
#include <map>
#include <set>
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
 * Checkpointing should be set up by including Global.h. Some more
 * checkpointing functions (for handling containers) are available by including
 * checkpoint_containers.h. */
namespace OM {
namespace interventions{
    struct ComponentId;
}
namespace util {
namespace checkpoint {
    
    const long DEFAULT_MAX_LENGTH = 2000;
    
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
    void validateListSize (long length, long max = DEFAULT_MAX_LENGTH);
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
    
} } }	// end of namespaces
#endif
