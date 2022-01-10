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
    DecayFuncHet(double tMult): tMult(tMult) {}
public:
    /** Default value: should make all eval() calls return 0 (i.e. infinitely
     * old deployment). */
    DecayFuncHet() : tMult( numeric_limits<double>::infinity() ) {}
    
    inline double getTMult() const{
        return tMult;
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        tMult & stream;
    }
    
    friend class BaseHetDecayFunction;
    friend class ConstantDecayFunction;
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
    static unique_ptr<DecayFunction> makeObject(
        const scnXml::DecayFunction& elt, const char* eltName
    );
    /** Return an object representing no decay (useful default). */
    static unique_ptr<DecayFunction> makeConstantObject();
    
    /** Return a value in the range [0,1] describing remaining effectiveness of
     * the intervention.
     * 
     * @param age Age of intervention/decayed property
     * @param sample A DecayFuncHet value sampled for the intervention and
     *  individual.
     * 
     * NOTE: As it is, values are calculated for the end of the time-period
     * being updated over. It would be more accurate to return the mean value
     * over this period (from age-1 to age), but difference should be small for
     * interventions being effective for a month or more. */
    inline double eval( SimTime age, DecayFuncHet sample )const{
        return eval( age * sample.getTMult() );
    }
    
    /** Sample a DecayFuncHet value (should be stored per individual).
     * 
     * Note that a DecayFuncHet is needed to call eval() even if heterogeneity
     * is not wanted. If sigma = 0 then the random number stream will not be
     * touched. */
    virtual DecayFuncHet hetSample (LocalRng& rng) const =0;
    
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
    virtual SimTime sampleAgeOfDecay (LocalRng& rng) const =0;
    
protected:
    DecayFunction() {}
    // Protected version. Note that the het sample parameter is needed even
    // when heterogeneity is not used so don't try calling this without that.
    virtual double eval(double ageDays) const =0;
};

} }
#endif
