/*
 This file is part of OpenMalaria.

 Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

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

#ifndef Hmod_DecayFunction
#define Hmod_DecayFunction

#include "Global.h"
#include <limits>
#include <boost/shared_ptr.hpp>

namespace scnXml
{
    class DecayFunction;
    class DecayFunctionValue;
}

using boost::shared_ptr;

namespace OM {
namespace util {

/** An interface for a few types of decay function (some of which may also be
 * suitible survival functions).
 *
 * Implemented functions:
 *
 * No heterogeneity implemented here; no need for an instance per individual.
 *****************************************************************************/
class DecayFunction
{
public:
    virtual ~DecayFunction() {}
    
    /** Return a new decay function, constructed from an XML element.
     * 
     * @param elt XML element specifying which function to use and parameters
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    static shared_ptr<DecayFunction> makeObject(
        const scnXml::DecayFunction& elt, const char* eltName
    );
    /** Return an object representing no decay (useful default). */
    static shared_ptr<DecayFunction> makeConstantObject();
    
    /** Return a value for age ageTS in time steps. */
    virtual double eval(TimeStep age) const =0;
    
protected:
    DecayFunction() {}
};

/** Wrapper around DecayFunction adding initial value. */
class DecayFunctionValue
{
    double initial;
    shared_ptr<DecayFunction> decayFunc;
public:
    DecayFunctionValue() : initial(numeric_limits<double>::quiet_NaN()) {}
    
    /** Assignment from XML element. */
    void set (const scnXml::DecayFunctionValue& elt, const char* eltName);
    
    /** Return true if decay function was never initialized. */
    bool notSet (){
        return decayFunc.get() == 0;
    }
    
    /** Return value for age ageTS in time steps. */
    double eval(TimeStep age) const{
        return initial * decayFunc->eval(age);
    }
};

} }
#endif