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
#ifndef OM_INTERVENTIONS_COHORT
#define OM_INTERVENTIONS_COHORT

// The includes here are more for documentation than required.
#include "interventions/HumanComponents.h"
#include <schema/interventions.h>

namespace OM { namespace interventions {

class CohortSelectionEffect : public HumanInterventionEffect {
public:
    CohortSelectionEffect( EffectId, const scnXml::Cohort& cohort );
    virtual ~CohortSelectionEffect() {}
    
    virtual void deploy( Host::Human& human, Deployment::Method method, VaccineLimits ) const;
    
    virtual Effect::Type effectType() const;
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const;
#endif
    
    /** For each RemoveAtCode (excluding NUM), this is a list of all cohort Ids
     * for which the option is enabled. */
    static vector<EffectId> removeAtIds[Cohort::REMOVE_AT_NUM];
};

} }
#endif
