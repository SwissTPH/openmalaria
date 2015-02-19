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

#include "PkPd/Drug/LSTMDrugThreeComp.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <boost/math/constants/constants.hpp>

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrugThreeComp::LSTMDrugThreeComp(const LSTMDrugType& type) :
    LSTMDrug(type.sample_Vd()),
    typeData(type),
    concA(0.0), concB(0.0), concC(0.0) /*, concABC(0.0)*/
{
    // These are from the Monolix article, pp38-39.
    const double k = type.sample_k();
    const double k12 = type.sample_k12(), k21 = type.sample_k21();
    const double k13 = type.sample_k13(), k31 = type.sample_k31();
    const double a0 = k*k21*k31;
    const double a1 = k*k31 + k21*k31 + k21*k13 + k*k21 + k31*k12;
    const double a2 = k + k12 + k13 + k21 + k31;
    const double third = 1.0 / 3.0;
    const double p = a1 - third * a2 * a2;
    const double q = (2.0 / 27.0) * a2 * a2 * a2 - third * a1 * a2 + a0;
    const double at = third * a2;
    const double rt = -third * p;
    const double r2 = 2.0 * sqrt(rt);
    const double phi = acos(-q / (rt * r2)) * third;
    const double twopithree = 2.0 * third * boost::math::constants::pi<double>();
    alpha = at - r2 * cos(phi);
    beta = at - r2 * cos(phi + twopithree);
    gamma = at - r2 * cos(phi + 2.0 * twopithree);
    const double invV = 1.0 / type.sample_Vd();
    // this assumes instantaneous absorbtion: k_a → ∞
    A = invV * (k21 - alpha) * (k31 - alpha) / ((alpha - beta) * (alpha - gamma));
    B = invV * (k21 - beta) * (k31 - beta) / ((beta - alpha) * (beta - gamma));
    C = invV * (k21 - gamma) * (k31 - gamma) / ((gamma - beta) * (gamma - alpha));
}

size_t LSTMDrugThreeComp::getIndex() const {
    return typeData.getIndex();
}
double LSTMDrugThreeComp::getConcentration() const {
    //NOTE: assuming no ABC term (see declaration of concA, concB, etc.).
    return concA + concB + concC /*- concABC*/;
}

void LSTMDrugThreeComp::medicate(double time, double qty, double bodyMass)
{
    medicate_vd(time, qty, vol_dist * bodyMass);
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrugThreeComp::calculateDrugFactor(uint32_t genotype) const {
    if( getConcentration() == 0.0 && doses.size() == 0 ) return 1.0; // nothing to do
    
    //FIXME: calculate factor
    return 1.0;
}

void LSTMDrugThreeComp::updateConcentration () {
    if( getConcentration() == 0.0 && doses.size() == 0 ) return;     // nothing to do
    
    // exponential decay of existing quantities:
    //TODO: is it faster to pre-calculate these and either store extra
    // parameters or adapt uses of alpha, beta, gamma, etc. below?
    concA *= exp(-alpha);
    concB *= exp(-beta);
    concC *= exp(-gamma);
    
    size_t doses_taken = 0;
    typedef pair<double,double> TimeConc;
    foreach( TimeConc& time_conc, doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            // add dose (instantaneous absorbtion):
            concA += A * time_conc.second * exp(-alpha * (1.0 - time_conc.first));
            concB += B * time_conc.second * exp(-beta * (1.0 - time_conc.first));
            concC += C * time_conc.second * exp(-gamma * (1.0 - time_conc.first));
            doses_taken += 1;
        }else /*i.e. tomorrow or later*/{
            time_conc.first -= 1.0;
        }
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    doses.erase(doses.begin(), doses.begin() + doses_taken);
    
    util::streamValidate( getConcentration() );
    if( getConcentration() < typeData.getNegligibleConcentration() ){
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        concA = 0.0;
        concB = 0.0;
        concC = 0.0;
        /*concABC = 0.0;*/
    }
}

void LSTMDrugThreeComp::checkpoint(ostream& stream){
    concA & stream;
    concB & stream;
    concC & stream;
    /*concABC & stream;*/
}
void LSTMDrugThreeComp::checkpoint(istream& stream){
    concA & stream;
    concB & stream;
    concC & stream;
    /*concABC & stream;*/
}

}
}