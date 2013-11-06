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

#include "interventions/Cohort.h"
#include "Host/Human.h"

namespace OM { namespace interventions {

vector<size_t> CohortSelectionEffect::removeAtIds[CohortSelectionEffect::REMOVE_AT_NUM];

CohortSelectionEffect::CohortSelectionEffect( size_t index, const scnXml::Cohort& cohort ) :
        HumanInterventionEffect(index)
{
    if( cohort.getFirstBoutOnly() ){
        removeAtIds[REMOVE_AT_FIRST_BOUT].push_back( index );
    }
    if( cohort.getFirstInfectionOnly() ){
        removeAtIds[REMOVE_AT_FIRST_INFECTION].push_back( index );
    }
    if( cohort.getFirstTreatmentOnly() ){
        removeAtIds[REMOVE_AT_FIRST_TREATMENT].push_back( index );
    }
}

void CohortSelectionEffect::deploy( Host::Human& human, Deployment::Method method )const{
    human.addToCohort( getIndex() );
}

Effect::Type CohortSelectionEffect::effectType() const{
    return Effect::COHORT;
}

} }
