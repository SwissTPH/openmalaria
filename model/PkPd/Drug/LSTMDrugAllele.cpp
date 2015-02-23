/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

// Declaration in LSTMDrugType.h due to circular dependency:
#include "PkPd/Drug/LSTMDrugType.h"
#include "util/errors.h"
#include "schema/pharmacology.h"

#include <boost/functional/hash.hpp>

using namespace std;

namespace OM { namespace PkPd {

//TODO: I think this Cache thing can be removed, though it may be useful as a template for 3-compartment/conversion models
LSTMDrugPD::Cache::Cache( double c, double d, double r ) :
    C0(c), duration(d), rate(r),
    C1(numeric_limits<double>::signaling_NaN()),
    drugFactor(numeric_limits<double>::signaling_NaN())
{
    // Generate hash using XOR and boost::hash
    boost::hash<double> hasher;
    hash = hasher(c) ^ hasher(d) ^ hasher(r);
}

LSTMDrugPD::LSTMDrugPD( const scnXml::Phenotype& phenotype ){
    n = phenotype.getSlope ();
    Kn = pow(phenotype.getIC50 (), n);
    V = phenotype.getMax_killing_rate ();  
    if( phenotype.getIC50().getSigma() > 0.0 ){
        throw util::unimplemented_exception("sampling IC50");
    }
}

double LSTMDrugPD::calcFactor( double neg_elim_rate, double& C0, double duration ) const{
    const double C1 = C0 * exp(neg_elim_rate * duration);
    
    // From Hastings & Winter 2011 paper
    // Note: these look a little different from original equations because Kn
    // is calculated when parameters are read from the scenario document instead of now.
    const double numerator = Kn + pow(C1, n);
    const double denominator = Kn + pow(C0, n);
    
    //TODO(performance): can we cache the value for each parameter combination?
    C0 = C1;    // C0 is an in/out parameter
    const double power = V / (-neg_elim_rate * n);
    return pow( numerator / denominator, power );       // unitless
}

} }
