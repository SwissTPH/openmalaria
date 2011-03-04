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
class DecayFunctionSuite;

using boost::shared_ptr;

namespace OM {
namespace util {

/** A sample of parameters used to make decay functions heterogenious.
 *
 * The default constructor only sets an NaN value; new instances must be
 * sampled by DecayFunction::hetSample() before use. */
class DecayFuncHet {
    double tMult;
public:
    DecayFuncHet() : tMult( numeric_limits<double>::quiet_NaN() ) {}
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        tMult & stream;
    }
    
    friend class BaseHetDecayFunction;
    friend class StepDecayFunction;
    friend class LinearDecayFunction;
    friend class ExponentialDecayFunction;
    friend class WeibullDecayFunction;
    friend class HillDecayFunction;
    friend class SmoothCompactDecayFunction;
    friend class ::DecayFunctionSuite;
 };

/** An interface for a few types of decay function (some of which may also be
 * suitible survival functions).
 *
 * Heterogeneity is implemented by passing a DecayFuncHet object to the eval
 * function; this should be sampled by the same DecayFunction.
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
    
    /** Sample a DecayFuncHet value (should be stored per individual). */
    virtual DecayFuncHet hetSample () const =0;
    
    /** Return a value in the range [0,1] describing remaining effectiveness of
     * the intervention.
     * 
     * @param age Age of intervention/decayed property.
     * @param sample A DecayFuncHet value sampled for the intervention and
     *  individual. */
    virtual double eval(TimeStep age, DecayFuncHet sample) const =0;
    
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
    
    /** As DecayFunction::eval(), but multiplied by initial value. */
    double eval(TimeStep age, DecayFuncHet sample) const{
        return initial * decayFunc->eval(age, sample);
    }
};

} }
#endif