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

class CohortSelectionComponent : public HumanInterventionComponent {
public:
    CohortSelectionComponent( ComponentId, const scnXml::Cohort& cohort );
    virtual ~CohortSelectionComponent() {}
    
    virtual void deploy( Host::Human& human, Deployment::Method method, VaccineLimits ) const;
    
    virtual Component::Type componentType() const;
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const;
#endif
    
    /** Return whether or not any of the map of component identifiersr passed
     * is a cohort identifier.
     * 
     * The map passed is just used like a set; its values are ignored. */
    static bool inAnyCohort( const map<ComponentId,TimeStep>& components );
};

} }
#endif
