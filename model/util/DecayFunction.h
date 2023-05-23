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

namespace OM {
namespace util {

/** An interface for a few types of decay function (some of which may also be
 * suitible survival functions).
 *
 * Heterogeneity is implemented by passing a DecayFunctionHet object to the eval
 * function; this should be sampled by the same DecayFunction.
 *****************************************************************************/
class DecayFunction
{
public:
    DecayFunction(bool increasing = false) : increasing(increasing) {}

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
    virtual unique_ptr<DecayFunction> hetSample (LocalRng& rng) const =0;
    
    /** Generate a DecayFunctionHet value from an existing sample. */
    virtual unique_ptr<DecayFunction> hetSample (NormalSample sample) const =0;
    
    /** Say you have a population of objects which each have two states:
     * decayed and not decayed. If you want to use a DecayFunction to model
     * the proportion of objects which have decayed, you need to work out per
     * object the age of decay. This function does that.
     * 
     * This is only valid where mu and sigma parameters are zero.
     * 
     * @returns Age at which an object should decay. */
    virtual SimTime sampleAgeOfDecay (LocalRng& rng) const =0;
    
    inline double eval(double ageDays) const {
        if(increasing)
            return 1.0 - _eval( ageDays );
        else
            return _eval( ageDays );
    }

        /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        // timeFactorHet & stream;
    }

protected:
    virtual double _eval(double ageDays) const =0;

    bool increasing;
};

} }
#endif
