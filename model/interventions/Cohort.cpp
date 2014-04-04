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

#include "interventions/Cohort.h"
#include "Host/Human.h"
#include "Monitoring/Survey.h"

namespace OM { namespace interventions {
    using namespace Monitoring::Report;

std::set<ComponentId> CohortSelectionComponent::cohortComponents;

CohortSelectionComponent::CohortSelectionComponent( ComponentId id, const scnXml::Cohort& cohort ) :
        HumanInterventionComponent(id, MI_NUM_ADDED_COHORT, MI_NUM_ADDED_COHORT)
{
    cohortComponents.insert( id );
}

void CohortSelectionComponent::deploy( Host::Human& human, Deployment::Method method, VaccineLimits )const{
    // Note: from the point view of the cohort, it makes sense to reset
    // healthSystemMemory. This was previously done by calling flushReports(),
    // but this resets health system memory for all reports. Instead we do not
    // reset, which is considered an acceptable approximation.
    
    //TODO(monitoring): reporting is inappropriate
    Monitoring::Survey::current().addInt( reportMeasure(method), human, 1 );
}

Component::Type CohortSelectionComponent::componentType() const{
    return Component::COHORT;
}
    
#ifdef WITHOUT_BOINC
void CohortSelectionComponent::print_details( std::ostream& out )const{
    out << id().id << "\tcohort selection";
}
#endif

} }
