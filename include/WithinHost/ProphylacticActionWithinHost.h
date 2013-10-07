/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2013 Swiss Tropical and Public Health Institute 
 * Copyright (C) 2005-2013 Liverpool School Of Tropical Medicine
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

#ifndef OM_WITHINHOST_PROPHYLACTIC_ACTION
#define OM_WITHINHOST_PROPHYLACTIC_ACTION

#include "WithinHost/DescriptiveWithinHost.h"

namespace OM { namespace WithinHost {
    
/** Extension to the DescriptiveWithinHostModel, including prophylactic
 * action of drugs. Partial alternative to the IPT model.
 *
 * NOTE: This shouldn't be implemented as a withinhost model, it should be a
 * separate model. For now this is easiest since it mirrors how the IPT code
 * works. */
class ProphylacticActionWithinHost : public DescriptiveWithinHostModel {
public:
    ProphylacticActionWithinHost () {}
    
    virtual void addProphylacticEffects(const vector<double>& pClearanceByTime);
    
protected:
    virtual void drugAction();
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
private:
    typedef list<pair<double,uint32_t> > Pending_t;
    /// A list of clearance chances by day, along with the number of factors acting on that day.
    Pending_t pendingClearanceProbabilities;
};

} }
#endif