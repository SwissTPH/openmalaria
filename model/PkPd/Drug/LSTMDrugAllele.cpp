/*
  This file is part of OpenMalaria.
 
  Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
  OpenMalaria is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or (at
  your option) any later version.
 
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

// Declaration in LSTMDrugType.h due to circular dependency:
#include "PkPd/Drug/LSTMDrugType.h"
#include "inputData.h"

#include <cmath>
#include <stdexcept>
#include <boost/functional/hash.hpp>

#include <gsl/gsl_integration.h>

using namespace std;

namespace OM { namespace PkPd {

LSTMDrugAllele::Cache::Cache( double c, double d, double r ) :
    C0(c), duration(d), rate(r) {
    assert(sizeof(double)==sizeof(hash));
    // Generate hash using XOR and boost::hash
    boost::hash<double> hasher;
    hash = hasher(c) ^ hasher(d) ^ hasher(r);
}

LSTMDrugAllele::LSTMDrugAllele( const scnXml::Allele& allele, double elimination_rate_constant ){
    slope = allele.getSlope ();
    power = allele.getMax_killing_rate () / (elimination_rate_constant * slope);
    IC50_pow_slope = pow(allele.getIC50 (), slope);
    max_killing_rate = allele.getMax_killing_rate ();  
}

double LSTMDrugAllele::calcFactor( const LSTMDrugType& drug, double& C1, double duration ) const{
    double C0 = C1;
    drug.updateConcentration( C1, duration );
    
    // Note: these look a little different from original equations because IC50_pow_slope
    // and power are calculated when read from the scenario document instead of here.
    double numerator = IC50_pow_slope + pow(C1, slope);
    double denominator = IC50_pow_slope + pow(C0, slope);
    
    return pow( numerator / denominator, power );
}

struct IV_conc_params {
    double C0;
    double ivRate; // dose->second.qty
    double neg_elimination_rate_constant;
    double elim_rate_dist;      // -neg_elimination_rate_constant * vol_dist
    double slope;
    double max_kill_rate;
    double IC50_pow_slope;
};
/** Function for calculating killing factor of IV and requiring integration.
 */
double func_IV_conc( double t, void* pp ){
    IV_conc_params *p = static_cast<IV_conc_params*>( pp );
    
    double conc_decay = exp(p->neg_elimination_rate_constant * t);
    double infusion = p->ivRate * (1.0 - conc_decay ) / p->elim_rate_dist
        + p->C0 * conc_decay;
    double infusion_pow_slope = pow(infusion, p->slope);
    double fC = p->max_kill_rate * infusion_pow_slope / (infusion_pow_slope + p->IC50_pow_slope);
    return fC;
}

double LSTMDrugAllele::calcFactorIV( const LSTMDrugType& drug, double& C0, double duration, double rate ) const{
    //TODO: properly test cache function and performance impact once we have a test scenario
    Cache key( C0, duration, rate );
    CachedIV::const_iterator it = cachedIV.find( key );
    if( it != cachedIV.end() ){
        // NOTE: we assume hash-conflicts won't occur, although this is not impossible
        // So do a little check (TODO: maybe make this a debug-only assert later?)
        if( !(key == *it) ){
            throw runtime_error( "LSTMDrugAllele::calcFactorIV: hash collision" );
        }
        
        C0 = it->C1;
        return it->drugFactor;
    } else {
        IV_conc_params p;
        p.C0 = C0;
        p.ivRate = rate;
        p.neg_elimination_rate_constant = drug.neg_elimination_rate_constant;
        p.elim_rate_dist =-p.neg_elimination_rate_constant * drug.vol_dist;
        p.slope = slope;
        p.max_kill_rate = max_killing_rate;
        p.IC50_pow_slope = IC50_pow_slope;
        
        gsl_function F;
        F.function = &func_IV_conc;
        F.params = static_cast<void*>(&p);
        
        //TODO: consider desired limits
        double abs_eps = 0.0, rel_eps = 1e-7;
        double intfC, err_eps;
        size_t n_evals;
        
        gsl_integration_qng (&F, 0.0, duration, abs_eps, rel_eps, &intfC, &err_eps, &n_evals); 
        
#ifndef NDEBUG
        cout << "IV: iterations: "<<n_evals<<"; error epsilon: "<<err_eps<<endl;
#endif
        //TODO: check err_eps is OK
        //TODO: check qng is the best integration algorithm
        
        drug.updateConcentrationIV( C0, duration, rate );
        
        // use key as our new cache
        key.C1 = C0;
        key.drugFactor = 1.0 / exp( intfC );
        bool inserted = cachedIV.insert( key ).second;
        assert( inserted );
        
        return key.drugFactor;
    }
}

} }
