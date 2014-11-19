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

#ifndef Hmod_Survey
#define Hmod_Survey

#include "Global.h"
#include "WithinHost/Diagnostic.h"
#include "Parameters.h"
#include "util/checkpoint_containers.h"
#include <bitset>
#include <map>

class UnittestUtil;
namespace scnXml{ class Monitoring; }
namespace OM {
namespace Host {
    class Human;
}
namespace Monitoring {
using WithinHost::Diagnostic;

/// Data struct for a single survey.
class Survey {
public:
    /// Initialize static parameters.
    static void init( const OM::Parameters& parameters,
                   const scnXml::Scenario& scenario,
                   const scnXml::Monitoring& monitoring );
    
    ///@brief Static access functions: these are here so that most users don't need to include Surveys.h
    //@{
    static inline const Diagnostic& diagnostic(){
        assert( m_diagnostic != 0 );
        return *m_diagnostic;
    }
    
    /** Humans should store a "cohort set" identifier which is initially 0.
     * Whenever a human gains or loses membership status in some
     * sup-population, it should update that value with this function.
     * 
     * @param old       Old identifier value (initially 0)
     * @param subPop    Sub-population to which membership status changed
     * @param isMember  New membership status
     * @returns         New identifier value
     */
    static uint32_t updateCohortSet( uint32_t old,
        interventions::ComponentId subPop, bool isMember );
    //@}
    
private:
    // ———  static members  ———
    
    static const Diagnostic* m_diagnostic;
    
    friend class SurveysType;
    friend class ::UnittestUtil;
};

} }
#endif
