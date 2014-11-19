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

#include "Monitoring/Surveys.h"
#include "Host/Human.h"
#include "util/ModelOptions.h"
#include "util/errors.h"
#include "util/CommandLine.h"
#include "schema/monitoring.h"
#include "schema/scenario.h"    // TODO: only for analysisNo

#include <boost/static_assert.hpp>
#include <stdexcept>

namespace OM { namespace Monitoring {
using WithinHost::diagnostics;

// -----  Static members  -----

const Diagnostic* Survey::m_diagnostic = 0;


void Survey::init( const OM::Parameters& parameters,
                   const scnXml::Scenario& scenario,
                   const scnXml::Monitoring& monitoring ){
    const scnXml::Surveys& surveys = monitoring.getSurveys();
    if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
        // So far the implemented Vivax code does not produce parasite
        // densities, thus this diagnostic model cannot be used.
        m_diagnostic = &diagnostics::make_deterministic( numeric_limits<double>::quiet_NaN() );
    }else if( surveys.getDetectionLimit().present() ){
        if( surveys.getDiagnostic().present() ){
            throw util::xml_scenario_error( "monitoring/surveys: do not "
                "specify both detectionLimit and diagnostic" );
        }
        if( util::CommandLine::option( util::CommandLine::DEPRECATION_WARNINGS ) ){
            std::cerr << "Deprecation warning: monitoring/surveys: "
                "specification of \"diagnostic\" is suggested over \"detectionLimit\"" << std::endl;
        }
        
        // This controls whether the detection limit is specified relative to
        // the Garki or other methods.
        double densitybias = numeric_limits<double>::quiet_NaN();
        if (util::ModelOptions::option (util::GARKI_DENSITY_BIAS)) {
            densitybias = parameters[Parameters::DENSITY_BIAS_GARKI];
        } else {
            if( scenario.getAnalysisNo().present() ){
                int analysisNo = scenario.getAnalysisNo().get();
                if ((analysisNo >= 22) && (analysisNo <= 30)) {
                    cerr << "Warning: these analysis numbers used to mean "
                        "use Garki density bias. If you do want to use this, "
                        "specify the option GARKI_DENSITY_BIAS; if not, nothing's wrong." << endl;
                }
            }
            densitybias = parameters[Parameters::DENSITY_BIAS_NON_GARKI];
        }
        double detectionLimit = surveys.getDetectionLimit().get() * densitybias;
        m_diagnostic = &diagnostics::make_deterministic( detectionLimit );
    }else{
        if( !surveys.getDiagnostic().present() ){
            throw util::xml_scenario_error( "monitoring/surveys: require "
                "either detectionLimit or diagnostic" );
        }
        if( util::ModelOptions::option(util::GARKI_DENSITY_BIAS) ){
            throw util::xml_scenario_error( "Use of GARKI_DENSITY_BIAS is not "
                "appropriate when monitoring/surveys/diagnostic is used." );
        }
        m_diagnostic = &diagnostics::get( surveys.getDiagnostic().get() );
    }
}
} }