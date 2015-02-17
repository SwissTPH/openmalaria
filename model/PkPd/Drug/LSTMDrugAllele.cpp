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

#include <cmath>
#include <stdexcept>
#include <boost/functional/hash.hpp>

#include <gsl/gsl_integration.h>
#include <fstream>

using namespace std;

namespace OM { namespace PkPd {

LSTMDrugPD::Cache::Cache( double c, double d, double r ) :
    C0(c), duration(d), rate(r),
    C1(numeric_limits<double>::signaling_NaN()),
    drugFactor(numeric_limits<double>::signaling_NaN())
{
    // Generate hash using XOR and boost::hash
    boost::hash<double> hasher;
    hash = hasher(c) ^ hasher(d) ^ hasher(r);
}

LSTMDrugPD::LSTMDrugPD( const scnXml::Phenotype& phenotype, double elimination_rate_constant ){
    slope = phenotype.getSlope ();
    power = phenotype.getMax_killing_rate () / (elimination_rate_constant * slope);
    IC50_pow_slope = pow(phenotype.getIC50 (), slope);
    max_killing_rate = phenotype.getMax_killing_rate ();  
}

double LSTMDrugPD::calcFactor( const LSTMDrugType& drug, double& C0, double duration ) const{
    // exponential decay of drug concentration
    const double C1 = C0 * exp(drug.getNegElimintationRateConst() * duration);
    
    // From Hastings & Winter 2011 paper
    // Note: these look a little different from original equations because IC50_pow_slope
    // and power are calculated when read from the scenario document instead of here.
    double numerator = IC50_pow_slope + pow(C1, slope);
    double denominator = IC50_pow_slope + pow(C0, slope);
    
    //TODO(performance): can we cache the value for each parameter combination?
    C0 = C1;    // C0 is an in/out parameter
    return pow( numerator / denominator, power );       // unitless
}

} }
