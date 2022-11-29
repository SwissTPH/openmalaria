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

#include "Global.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "util/UnitParse.h"
#include "schema/scenario.h"
#include "mon/management.h"

namespace OM {

void sim::init( const scnXml::Scenario& scenario ){
    SimData::interval = scenario.getModel().getParameters().getInterval();
    SimData::steps_per_year = sim::inSteps(sim::oneYear());
    SimData::years_per_step = 1.0 / SimData::steps_per_year;
    sim::s_max_human_age = sim::fromYearsD( scenario.getDemography().getMaximumAgeYrs() );
    
    sim::s_start = sim::origin();
    auto mon = scenario.getMonitoring();
    if( mon.getStartDate().present() ){
        try{
            // on failure, this throws or returns sim::never()
            sim::s_start = UnitParse::parseDate( mon.getStartDate().get() );
            if( sim::s_start == sim::never() ){
                throw util::format_error( "invalid format (expected YYYY-MM-DD)" );
            }
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("monitoring/startDate: ").append(e.message()) );
        }
    }
    
    sim::s_interv = sim::never(); // large negative number
    sim::s_end = mon::readSurveyDates( mon );

    cout << "END: " << sim::s_end << endl;
}

}
