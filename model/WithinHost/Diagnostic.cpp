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
#include "util/random.h"
#include "util/errors.h"

namespace OM { namespace WithinHost {

Diagnostic Diagnostic::default_;

void Diagnostic::setDeterministic(double limit){
    assert( (boost::math::isnan)(dens_lim) );       // multiple initialisations
    specificity = numeric_limits<double>::quiet_NaN();
    dens_lim = limit;
}

void Diagnostic::setXml( const scnXml::HSDiagnostic& elt ){
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
                    string("diagnostic: specificity must be in range [0,1]") );
            }
        }
    }else{
        // This should be impossible since according to schema one of these
        // elements must be present.
        throw SWITCH_DEFAULT_EXCEPTION;
    }
    if( dens_lim < 0.0 ){
        throw util::xml_scenario_error(
            string("diagnostic: must have density ≥ 0") );
    }
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

} }
