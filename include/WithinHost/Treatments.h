/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_WithinHost_Treatments
#define Hmod_WithinHost_Treatments

#include "WithinHost/WHInterface.h"

#include <boost/ptr_container/ptr_vector.hpp>

using namespace std;

namespace scnXml{
    class TreatmentOption;
}

namespace OM {
namespace WithinHost {

/**
 * Objects of this class describe effects of a treatment, after selecting a
 * compliance/adherence/... group.
 * 
 * The static part provides a way to find objects.
 * 
 * For use within WithinHost only.
 */
class Treatments {
public:
    /// @brief Static methods
    //@{
    /// Add the "legacy" option
    static void init();
    
    /** Configure a new treatment option, and return the code used to select
     * that option later. */
    static TreatmentId addTreatment( const scnXml::TreatmentDescription& desc );
    
    /** Find an appropriate TreatmentGroup. If multiple groups exist for the
     * given TreatmentId, this will randomly select one. */
    static const Treatments& select( TreatmentId treatId );
    //@}
    
    /// @brief Types used by the non-static methods
    //@{
    /// Stages effected
    enum Stages{
        NONE /*i.e. no effect*/, LIVER, BLOOD, BOTH
    };
    /// Details of effects
    struct Action{
        Action(int len, Stages s): timesteps(len), stage(s) {}
        util::TimeStep timesteps;
        Stages stage;
    };
    //@}
    
    /// @brief Non-static methods
    //@{
    /** Get the list of actions.
     * 
     * There can be any number of actions; most likely are zero, one or two. */
    inline const vector<Action>& getEffects() const{ return effects; }
    //@}
    
private:
    // static:
    static boost::ptr_vector<std::vector<Treatments> > groups;
    
    // non-static:
    Treatments( double pCum );
    Treatments( double pCum, const scnXml::TreatmentOption& elt );
    
    double pCum;        // cumulative probability of selection; last should be 1
    // we can have several effects; each applies independently but simultaneously
    vector<Action> effects;
    
    friend class ::UnittestUtil;
};

}
}
#endif
