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

#ifndef Hmod_SimpleDecayingValue
#define Hmod_SimpleDecayingValue

#include "Global.h"
#include "util/DecayFunction.h"
#include <boost/shared_ptr.hpp>

namespace OM {
namespace util {

/** Parameters to for a decaying intervention whose effect can be described
 * by a single scalar output.
 * 
 * This class outputs a scalar which defaults to and decays to 0, but
 * immediately after deployment will be the set _initial value_.
 *****************************************************************************/
class SimpleDecayingValue
{
public:
    /** Default construction: always return 0. */
    SimpleDecayingValue() :
            initial(0.0),
            deploy_t(SimTime::never()) {}
    
    /** Configure from an XML element. */
    inline void set (double initial_value, const scnXml::DecayFunction& elt, const char* name){
        decay = DecayFunction::makeObject(elt, name);
        initial = initial_value;
    }
    
    /** Trigger a deployment, if decay function was set. */
    inline void deploy (SimTime time){
        if( decay.get() == 0 ) return;  // cannot deploy
        deploy_t = time;
        het = decay->hetSample();
    }
    
    /** Get the value (0 if before any deployment or after complete decay,
     * also 0 if no decay function or initial value was set,
     * otherwise between zero and the inital value). */
    inline double current_value (SimTime time) const{
        if( decay.get() == 0 ) return 0.0;  // decay wasn't set: the always return 0
        
        return initial * decay->eval( time - deploy_t, het );
    }
    
    /** Checkpointing: only checkpoint parameters which change after initial
     * set-up. */
    template<class S>
    void operator& (S& stream) {
        het & stream;
        deploy_t & stream;
    }
    
private:
    /** Description of decay of effects on emergence. */
    //C++11: could use unique_ptr and implement move support
    // For C++98 we cannot use auto_ptr because the class object needs to be stored in a vector
    boost::shared_ptr<util::DecayFunction> decay;
    
    /** Initial value. Is initialised to 0. */
    double initial;
    
    /** Description of larviciding decay. */
    util::DecayFuncHet het;
    
    /** Time of larviciding deployment. */
    SimTime deploy_t;
};

} }
#endif
