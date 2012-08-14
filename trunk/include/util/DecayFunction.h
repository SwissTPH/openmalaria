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
#include "util/sampler.h"
#include <limits>
#include <memory>

namespace scnXml
{
    class DecayFunction;
    class DecayFunctionValue;
}
class DecayFunctionSuite;

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
    
    inline double getTMult() const{
        return tMult;
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        tMult & stream;
    }
    
    friend class BaseHetDecayFunction;
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
    static auto_ptr<DecayFunction> makeObject(
        const scnXml::DecayFunction& elt, const char* eltName
    );
    /** Return an object representing no decay (useful default). */
    static auto_ptr<DecayFunction> makeConstantObject();
    
    /** Return a value in the range [0,1] describing remaining effectiveness of
     * the intervention.
     * 
     * @param age Age of intervention/decayed property during the time-step
     *  for which values are now being calculated.
     * @param sample A DecayFuncHet value sampled for the intervention and
     *  individual.
     * 
     * NOTE: As it is, values are calculated for the end of the time-period
     * being updated over. It would be more accurate to return the mean value
     * over this period (from age-1 to age), but difference should be small for
     * interventions being effective for a month or more. */
    virtual double eval(TimeStep age, DecayFuncHet sample) const =0;
    
    /** Sample a DecayFuncHet value (should be stored per individual). */
    virtual DecayFuncHet hetSample () const =0;
    
    /** Generate a DecayFuncHet value from an existing sample. */
    virtual DecayFuncHet hetSample (NormalSample sample) const =0;
    
    /** Say you have a population of objects which each have two states:
     * decayed and not decayed. If you want to use a DecayFunction to model
     * the proportion of objects which have decayed, you need to work out per
     * object the age of decay. This function does that.
     * 
     * This is only valid where mu and sigma parameters are zero.
     * 
     * @returns Age at which an object should decay. */
    virtual TimeStep sampleAgeOfDecay () const =0;
    
protected:
    DecayFunction() {}
};

} }
#endif