/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
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

class DecayFunctionHet;

/** An interface for a few types of decay function (some of which may also be
 * suitible survival functions).
 *
 * Heterogeneity is implemented by passing a DecayFunctionHet object to the eval
 * function; this should be sampled by the same DecayFunction.
 *****************************************************************************/
class DecayFunction
{
public:
    DecayFunction(const scnXml::DecayFunction& elt) :
        complement(elt.getComplement())
    {}

    virtual ~DecayFunction() {}
    
    /** Return a new decay function, constructed from an XML element.
     * 
     * @param elt XML element specifying which function to use and parameters
     * @param eltName Name of XML element (for reasonable error reporting)
     */
    static unique_ptr<DecayFunction> makeObject(
        const scnXml::DecayFunction& elt, const char* eltName
    );

    /** Sample a DecayFunctionHet value (should be stored per individual).
     * 
     * Note that a DecayFunctionHet is needed to call eval() even if heterogeneity
     * is not wanted. If sigma = 0 then the random number stream will not be
     * touched. */
    virtual DecayFunctionHet hetSample (LocalRng& rng) const =0;
    
    /** Generate a DecayFunctionHet value from an existing sample. */
    virtual DecayFunctionHet hetSample (NormalSample sample) const =0;
    
    /** Say you have a population of objects which each have two states:
     * decayed and not decayed. If you want to use a DecayFunction to model
     * the proportion of objects which have decayed, you need to work out per
     * object the age of decay. This function does that.
     * 
     * This is only valid where mu and sigma parameters are zero.
     * 
     * @returns Age at which an object should decay. */
    virtual SimTime sampleAgeOfDecay (LocalRng& rng) const =0;
    
    virtual double eval(double ageDays) const =0;

    bool complement;

protected:
    DecayFunction() : complement(true) {}
    // Protected version. Note that the het sample parameter is needed even
    // when heterogeneity is not used so don't try calling this without that.
};

/** A sample of parameters used to make decay functions heterogenious.
 *
 * The default constructor only sets an NaN value; new instances must be
 * sampled by DecayFunction::hetSample() before use. */
class DecayFunctionHet {
    double timeFactorHet;
public:
    DecayFunctionHet(double timeFactorHet, shared_ptr<DecayFunction> d, bool complement = true): timeFactorHet(timeFactorHet), sample(d) {}

    /** Default value: should make all eval() calls return 0 (i.e. infinitely
     * old deployment). */
    DecayFunctionHet();
    
    inline double getTimeFactorHet() const{
        return timeFactorHet;
    }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        timeFactorHet & stream;
    }

    /** Return a value in the range [0,1] describing remaining effectiveness of
     * the intervention.
     * 
     * @param age Age of intervention/decayed property
     * @param sample A DecayFunctionHet value sampled for the intervention and
     *  individual.
     * 
     * NOTE: As it is, values are calculated for the end of the time-period
     * being updated over. It would be more accurate to return the mean value
     * over this period (from age-1 to age), but difference should be small for
     * interventions being effective for a month or more. */
    inline double eval( SimTime age )const
    {
        if(sample->complement)
            return 1.0 - sample->eval( age * timeFactorHet);
        else
            return sample->eval( age * timeFactorHet);
    }
    
    // friend class BaseHetDecayFunction;
    // friend class ConstantDecayFunction;

    // TEMPORARY
    shared_ptr<DecayFunction> sample;
 };

class NullDecayFunction : public DecayFunction
{
public:
    virtual double eval(double ageDays) const { return 0.0; }

    virtual DecayFunctionHet hetSample (LocalRng& rng) const { return DecayFunctionHet(); };
    
    /** Generate a DecayFunctionHet value from an existing sample. */
    virtual DecayFunctionHet hetSample (NormalSample sample) const { return DecayFunctionHet(); };

    virtual SimTime sampleAgeOfDecay (LocalRng& rng) const { return 0; };

};


} }
#endif
