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

#include <boost/math/constants/constants.hpp>
#include <gsl/gsl_integration.h>

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrugThreeComp::LSTMDrugThreeComp(const LSTMDrugType& type) :
    LSTMDrug(type.sample_Vd()),
    typeData(type),
    concA(0.0), concB(0.0), concC(0.0), concABC(0.0),
    nka(-type.sample_ka())
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
    const double pi23 = (2.0 / 3.0) * boost::math::constants::pi<double>();
    // negative alpha, beta, gamma:
    na = r2 * cos(phi) - at;
    nb = r2 * cos(phi + pi23) - at;
    ng = r2 * cos(phi + 2.0 * pi23) -at;
    // A, B, C from Monolix 1.3.3 (p44):
    const double kaV = -nka / type.sample_Vd();
    A = kaV * (k21 + na) * (k31 + na) / ((na - nka) * (nb - na) * (ng - na));
    B = kaV * (k21 + nb) * (k31 + nb) / ((nb - nka) * (na - nb) * (ng - nb));
    C = kaV * (k21 + ng) * (k31 + ng) / ((ng - nka) * (nb - ng) * (na - ng));
}

size_t LSTMDrugThreeComp::getIndex() const {
    return typeData.getIndex();
}
double LSTMDrugThreeComp::getConcentration() const {
    //NOTE: assuming no ABC term (see declaration of concA, concB, etc.).
    return concA + concB + concC - concABC;
}

void LSTMDrugThreeComp::medicate(double time, double qty, double bodyMass)
{
    medicate_vd(time, qty, vol_dist * bodyMass);
}

/// Parameters for func_fC
struct Params_fC {
    double cA, cB, cC, cABC;    // concentration parameters
    double na, nb, ng, nka;     // decay parameters
    double n;       // slope: unitless
    double V;       // max killing rate: unitless
    double Kn;      // IC50^n: (mg/kg) ^ n
};
/** Function for calculating concentration and then killing function at time t
 * 
 * @param t The variable being integrated over (in this case, time since start
 *      of day or last dose, units days)
 * @param pp Pointer to a Params_fC struct
 * @return killing rate (unitless)
 */
double func_fC( double t, void* pp ){
    const Params_fC& p = *static_cast<const Params_fC*>( pp );
    
    // exponential decay of drug concentration:
    const double concA = p.cA * exp(p.na * t);
    const double concB = p.cB * exp(p.nb * t);
    const double concC = p.cC * exp(p.ng * t);
    const double concABC = p.cABC * exp(p.nka * t);
    const double conc = concA + concB + concC - concABC;      // mg/l
    
    const double cn = pow(conc, p.n);        // (mg/l) ^ n
    const double fC = p.V * cn / (cn + p.Kn);       // unitless
    return fC;
}
double LSTMDrugThreeComp::calculateFactor(const Params_fC& p, double duration) const{
    gsl_function F;
    F.function = &func_fC;
    // gsl_function doesn't accept const; we re-apply const later
    F.params = static_cast<void*>(const_cast<Params_fC*>(&p));
    
    // TODO: verify these error limits are sufficient.
    const double abs_eps = 1e-3, rel_eps = 1e-3;
    const int qag_rule = 2; // GSL_INTEG_GAUSS21 should be sufficient?
    double intfC, err_eps;
    
    const size_t max_iterations = 1000;     // 100 seems enough, but no harm in using a higher value
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc (max_iterations);
    int r = gsl_integration_qag (&F, 0.0, duration, abs_eps, rel_eps, max_iterations, qag_rule, workspace, &intfC, &err_eps);
    if( r != 0 ){
        throw TRACED_EXCEPTION( "calcFactorIV: error from gsl_integration_qag",util::Error::GSL );
    }
    if( err_eps > 1e-2 ){
        // This could be a warning, except that warnings tend to be ignored.
        ostringstream msg;
        msg << "calcFactorIV: error epsilon is large: "<<err_eps<<" (integral is "<<intfC<<")";
        throw TRACED_EXCEPTION( msg.str(), util::Error::GSL );
    }
    gsl_integration_workspace_free (workspace); //TODO: don't always reallocate!
    
    return 1.0 / exp( intfC );  // drug factor
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrugThreeComp::calculateDrugFactor(uint32_t genotype) const {
    if( getConcentration() == 0.0 && doses.size() == 0 ) return 1.0; // nothing to do
    
    Params_fC p;
    p.cA = concA;       p.cB = concB;   p.cC = concC;   p.cABC = concABC;
    p.na = na;  p.nb = nb;      p.ng = ng;      p.nka = nka;
    const LSTMDrugPD& pd = typeData.getPD(genotype);
    p.n = pd.slope();   p.V = pd.max_killing_rate();    p.Kn = pd.IC50_pow_slope();
    
    double time = 0.0;  // time since start of day
    double totalFactor = 1.0;   // survival factor for whole day
    
    typedef pair<double,double> TimeConc;
    foreach( const TimeConc& time_conc, doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( time < time_conc.first ){
                double duration = time_conc.first - time;
                totalFactor *= calculateFactor(p, duration);
                p.cA *= exp(p.na * duration);
                p.cB *= exp(p.nb * duration);
                p.cC *= exp(p.ng * duration);
                p.cABC *= exp(p.nka * duration);
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add dose (instantaneous absorbtion):
            p.cA += A * time_conc.second;
            p.cB += B * time_conc.second;
            p.cC += C * time_conc.second;
            p.cABC += (A + B + C) * time_conc.second;
            // add dose (instantaneous absorbtion):
        }else /*i.e. tomorrow or later*/{
            // ignore
        }
    }
    if( time < 1.0 ){
        totalFactor *= calculateFactor(p, 1.0 - time);
    }
    
    return totalFactor;
}

void LSTMDrugThreeComp::updateConcentration () {
    if( getConcentration() == 0.0 && doses.size() == 0 ) return;     // nothing to do
    
    // exponential decay of existing quantities:
    //TODO(performance): is it faster to pre-calculate these and either store extra
    // parameters or adapt uses of alpha, beta, gamma, etc. below?
    concA *= exp(na);
    concB *= exp(nb);
    concC *= exp(ng);
    concABC *= exp(nka);
    
    size_t doses_taken = 0;
    typedef pair<double,double> TimeConc;
    foreach( TimeConc& time_conc, doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            // add dose (instantaneous absorbtion):
            concA += A * time_conc.second * exp(na * (1.0 - time_conc.first));
            concB += B * time_conc.second * exp(nb * (1.0 - time_conc.first));
            concC += C * time_conc.second * exp(ng * (1.0 - time_conc.first));
            concABC += (A + B + C) * time_conc.second * exp(nka * (1.0 - time_conc.first));
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
        concABC = 0.0;
    }
}

void LSTMDrugThreeComp::checkpoint(ostream& stream){
    concA & stream;
    concB & stream;
    concC & stream;
    concABC & stream;
}
void LSTMDrugThreeComp::checkpoint(istream& stream){
    concA & stream;
    concB & stream;
    concC & stream;
    concABC & stream;
}

}
}