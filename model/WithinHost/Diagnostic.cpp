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

#include "WithinHost/Diagnostic.h"
#include "Parameters.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/ModelOptions.h"
#include "schema/scenario.h"
#include <limits>
#include <boost/ptr_container/ptr_map.hpp>

namespace OM { namespace WithinHost {

// ———  Diagnostic (non-static)  ———

Diagnostic::Diagnostic( const Parameters& parameters, const scnXml::Diagnostic& elt ){
    if( elt.getDeterministic().present() ){
        specificity = numeric_limits<double>::quiet_NaN();
        dens_lim = elt.getDeterministic().get().getMinDensity();
    }else if( elt.getStochastic().present() ){
        dens_lim = elt.getStochastic().get().getDens_50();
        if( dens_lim == 0.0 ){
            // The equation used for stochastic diagnostics breaks down when
            // dens=dens_lim=0 and for other cases the deterministic model is
            // the same when dens_lim=0.
            specificity = numeric_limits<double>::quiet_NaN();
        }else{
            specificity = elt.getStochastic().get().getSpecificity();
            if( specificity < 0.0 || specificity > 1.0 ){
                throw util::xml_scenario_error(
                    string("diagnostics/diagnostic(").append(elt.getName())
                    .append("): specificity must be in range [0,1]") );
            }
        }
    }else{
        // This should be impossible since according to schema one of these
        // elements must be present.
        throw SWITCH_DEFAULT_EXCEPTION;
    }
    if( dens_lim < 0.0 ){
        throw util::xml_scenario_error(
            string("diagnostics/diagnostic(").append(elt.getName())
            .append("): must have density ≥ 0") );
    }
    
    // We use a bias factor to adjust the "units" used to specify the density
    // of this diagnostic, since estimates from Garki and the standard
    // non-Garki sources are not equivalent to those from the Malariatherapy
    // data (which is used internally).
    if( !elt.getUnits().present() ){
        if( util::ModelOptions::option(util::GARKI_DENSITY_BIAS) ){
            // User must be explicit in this case, because presumably the Garki
            // bias is to be used for some diagnostics but likely not all
            // (e.g. neonatal mortality).
            throw util::xml_scenario_error( "diagnostics/diagnostic(*)/units: must specify this attribute when GARKI_DENSITY_BIAS is set" );
        }
        // otherwise we assume "Other"
        dens_lim *= parameters[Parameters::DENSITY_BIAS_NON_GARKI];
    }else if( elt.getUnits().get() == "Other" ){
        dens_lim *= parameters[Parameters::DENSITY_BIAS_NON_GARKI];
    }else if( elt.getUnits().get() == "Garki" ){
        dens_lim *= parameters[Parameters::DENSITY_BIAS_GARKI];
    }else{
        assert( elt.getUnits().get() == "Malariatherapy" );
        // in this case we don't need to use a bias factor
    }
}

Diagnostic::Diagnostic(double minDens){
    specificity = numeric_limits<double>::quiet_NaN();
    dens_lim = minDens;
}

bool Diagnostic::isPositive( double dens ) const {
    if( (boost::math::isnan)(specificity) ){
        // use deterministic test
        return dens >= dens_lim;
    }else{
        // dens_lim is dens_50 in this case
        double pPositive = 1.0 + specificity * (dens / (dens + dens_lim) - 1.0);
//         double pPositive = (dens + dens_lim - dens_lim * specificity) / (dens + dens_lim);       // equivalent
        return util::random::bernoulli( pPositive );
    }
}


// ———  diagnostics (static)  ———

typedef boost::ptr_map<string,Diagnostic> Diagnostic_set;
Diagnostic_set diagnostic_set;

void diagnostics::clear(){
    diagnostic_set.clear();
}

void diagnostics::init( const Parameters& parameters, const scnXml::Diagnostics& diagnostics ){
    foreach( const scnXml::Diagnostic& diagnostic, diagnostics.getDiagnostic() ){
        string name = diagnostic.getName();     // conversion fails without this extra line
        bool inserted = diagnostic_set.insert( name, new Diagnostic(parameters, diagnostic) ).second;
        if( !inserted ){
            throw util::xml_scenario_error( string("diagnostic with this name already set: ").append(diagnostic.getName()) );
        }
    }
}

const Diagnostic& diagnostics::get( const string& name ){
    Diagnostic_set::const_iterator it = diagnostic_set.find(name);
    if( it == diagnostic_set.end() ){
        throw util::xml_scenario_error( string("diagnostic not found: ").append(name) );
    }
    return *it->second;
}

const Diagnostic& diagnostics::make_deterministic(double minDens){
    string name = "";   // anything not matching an existing name is fine
    assert( diagnostic_set.count(name) == 0 );  // we shouldn't need to call make_deterministic more than once, so I think this will do
    Diagnostic *diagnostic = new Diagnostic( minDens );
    diagnostic_set.insert( name, diagnostic );
    return *diagnostic;
}

} }
