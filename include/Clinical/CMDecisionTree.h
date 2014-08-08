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

#ifndef Hmod_ESDecisionTree
#define Hmod_ESDecisionTree

#include "Global.h"
#include "WithinHost/Pathogenesis/State.h"
#include "Clinical/Episode.h"
#include <schema/healthSystem.h>

namespace OM {
namespace WithinHost {
    class WHInterface;
}
namespace Clinical {
using WithinHost::WHInterface;

using std::auto_ptr;

/** All data which needs to be passed to the decision tree evaluators. */
struct CMHostData {
    CMHostData (double aY, WHInterface& wH, Episode::State pS) :
        ageYears(aY), withinHost(wH), pgState(pS) {}
    //TODO: do we need age here?
    double ageYears;
    WHInterface& withinHost;
    Episode::State pgState;
};


/**
 * Decision tree node abstraction.
 * 
 * Sub-classes represent either a decision node (first/second line case, a
 * diagnostic with positive/negative outcome, a random decision) or an action.
 *****************************************************************************/
class CMDecisionTree {
public:
    virtual ~CMDecisionTree() {}
    
    /// Create a user-configured decision from an XML node.
    static auto_ptr<CMDecisionTree> create( const ::scnXml::DecisionTree& node );
    
    /** Run the decision tree. */
    virtual void exec( CMHostData hostData ) const =0;
};

} }
#endif