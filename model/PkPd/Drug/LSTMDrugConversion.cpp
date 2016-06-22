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

#include "PkPd/Drug/LSTMDrugConversion.h"
#include "util/errors.h"
#include "util/StreamValidator.h"
#include "util/vectors.h"

#include <gsl/gsl_integration.h>
#include <limits>

using namespace std;

namespace OM { namespace PkPd {
    
LSTMDrugConversion::LSTMDrugConversion(const LSTMDrugType& parent,
                                       const LSTMDrugType& metabolite) :
    LSTMDrug(parent.sample_Vd()),
    parentType(parent),
    metaboliteType(metabolite),
    qtyG(0.0), qtyP(0.0), qtyM(0.0),
    nka(-parent.sample_ka()),
    nkP_sample(-parent.sample_elim_rate()),
    nconv_sample(-parent.sample_conv_rate()),
    nkM_sample(-metabolite.sample_elim_rate()),
    vol_dist_metabolite(metabolite.sample_Vd())
{}

//LSTMDrugConversion::~LSTMDrugConversion(){}

size_t LSTMDrugConversion::getIndex() const {
    return parentType.getIndex();       // parent drug's index should work
}
double LSTMDrugConversion::getConcentration(size_t index) const {
    if( index == parentType.getIndex() ){
        return getParentConcentration();
    }else if( index == metaboliteType.getIndex() ){
        return getMetaboliteConcentration();
    }else return 0.0;
}
double LSTMDrugConversion::getParentConcentration() const {
    return qtyP / (vol_dist * last_bm);
}
double LSTMDrugConversion::getMetaboliteConcentration() const {
    return qtyM / (vol_dist_metabolite * last_bm);
}


void LSTMDrugConversion::medicate(double time, double qty, double bodyMass)
{
    // Note: we use second value of dose pair as quanity (mg) not concentration
    // TODO: this is correct in AS unittest (qty = 200 = 4 * 50), but why is qty mg/kg not mg?
    medicate_vd(time, qty, 1.0);
}

/// Parameters for func_convFactor and LSTMDrugConversion::calculateFactor
struct Params_convFactor {
    // Quantities of parent in gut, parent in circulation and metabolite in
    // circulation, in mg (A, B, C in paper):
    double qtyG, qtyP, qtyM;
    // nka = -x (negative of absorption rate)
    // nkM = -k
    // nl = -(y + z)
    double nka, nkM, nl;
    // terms involving x, y, z and the molecular weight ratio:
    double f, g, h, i, j;
    
    double invVdP, invVdM;     // 1.0 / (Vd * body mass) for parent and metabolite; units: 1/l
    double nP, nM;       // slope: unitless
    double VP, VM;       // max killing rate: unitless
    double KnP, KnM;      // IC50^n: (mg/kg) ^ n
};

double calculateParentQuantity( struct Params_convFactor p, double expAbsorb, double expPLoss ) {
    return p.f * p.qtyG * expAbsorb
        + (p.qtyP - p.f * p.qtyG) * expPLoss;
}

double calculateParentDrugFactor( struct Params_convFactor p, double expAbsorb, double expPLoss ) {
    const double qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
    const double cP = qtyP * p.invVdP;                  // concentrations; mg/l*/
    const double cnP = pow(cP, p.nP);                   // (mg/l) ^ n
    const double fCP = p.VP * cnP / (cnP + p.KnP);      // unitless
    return fCP;
}

double calculateMetaboliteQuantity(struct Params_convFactor p, double expAbsorb, double expPLoss, double t) {
    return p.g * p.qtyG * expAbsorb
        + (p.h * p.qtyG - p.i * p.qtyP) * expPLoss
        + (p.j * p.qtyG + p.i * p.qtyP + p.qtyM) * exp(p.nkM * t);
}

double calculateMetaboliteDrugFactor( struct Params_convFactor p, double expAbsorb, double expPLoss, double t ) {
    const double qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, t);
    const double cM = qtyM * p.invVdM;              // concentrations; mg/l
    const double cnM = pow(cM, p.nM);               // (mg/l) ^ n
    const double fCM = p.VM * cnM / (cnM + p.KnM);  // unitless
    return fCM;
}

/** Function for calculating concentration and then killing function at time t
 * 
 * @param t The variable being integrated over (in this case, time since start
 *      of day or last dose, units days)
 * @param pp Pointer to a Params_convFactor struct
 * @return killing rate (unitless)
 */
double func_convFactor( double t, void* pp ){
    const Params_convFactor& p = *static_cast<const Params_convFactor*>( pp );
    const double expAbsorb = exp(p.nka * t), expPLoss = exp(p.nl * t);
    const double fCP = calculateParentDrugFactor( p, expAbsorb, expPLoss );
    const double fCM = calculateMetaboliteDrugFactor( p, expAbsorb, expPLoss, t );
    // use the most effective killing factor (from area under the drug kill curve), which is the one with the bigger number
    return max(fCP,fCM);
}

const size_t GSL_INTG_CONV_MAX_ITER = 1000;     // 10 seems enough, but no harm in using a higher value
gsl_integration_workspace *gsl_intgr_conv_wksp = gsl_integration_workspace_alloc (GSL_INTG_CONV_MAX_ITER);
//NOTE: we "should" free, but mem-leaks at end of program aren't really important
// gsl_integration_workspace_free (gsl_intgr_conv_wksp);
double LSTMDrugConversion::calculateFactor(const Params_convFactor& p, double duration) const{
    gsl_function F;
    F.function = &func_convFactor;
    // gsl_function doesn't accept const; we re-apply const later
    F.params = static_cast<void*>(const_cast<Params_convFactor*>(&p));
    
    // NOTE: tolerances are arbitrary, but seem to be sufficient
    const double abs_eps = 1e-3, rel_eps = 1e-3;
    // NOTE: 1 through 6 are different algorithms of increasing complexity
    const int qag_rule = 1;     // alg 1 seems to be good enough
    double intfC, err_eps;      // intfC will carry our result; err_eps is a measure of accuracy of the result
    
    int r = gsl_integration_qag (&F, 0.0, duration, abs_eps, rel_eps,
                                 GSL_INTG_CONV_MAX_ITER, qag_rule, gsl_intgr_conv_wksp, &intfC, &err_eps);
    if( r != 0 ){
        throw TRACED_EXCEPTION( "calculateFactor: error from gsl_integration_qag",util::Error::GSL );
    }
    if( err_eps > 5e-2 ){
        // This could be a warning, except that warnings tend to be ignored.
        ostringstream msg;
        msg << "calculateFactor: error epsilon is large: "<<err_eps<<" (integral is "<<intfC<<")";
        throw TRACED_EXCEPTION( msg.str(), util::Error::GSL );
    }
    return exp( -intfC );  // drug factor
}

void LSTMDrugConversion::setConversionParameters(Params_convFactor& p,
        double body_mass) const
{
    p.qtyG = qtyG; p.qtyP = qtyP; p.qtyM = qtyM;
    
    // decay "constants" (dependent on body mass):
    const double nkP = nkP_sample * pow(body_mass, parentType.neg_m_exponent());      // -y
    const double nconv = nconv_sample * pow(body_mass, parentType.neg_m_exponent());  // -z
    p.nka = nka;        // -x
    p.nkM = nkM_sample * pow(body_mass, metaboliteType.neg_m_exponent());  // -k
    p.nl = nkP + nconv;    // -(y + z)
    
    p.f = nka / (p.nl - nka);     // x*A' /  (y+z-x)
    const double rz = parentType.molecular_weight_ratio() * nconv;
    p.g = rz * nka / ((nka - p.nl) * (nka - p.nkM));
    p.h = rz * nka / ((nka - p.nl) * (p.nkM - p.nl));
    p.i = rz / (p.nl - p.nkM);
    p.j = rz * nka / ((p.nkM - p.nl) * (p.nkM - nka));
    
    p.invVdP = 1.0 / (vol_dist * body_mass); p.invVdM = 1.0 / (vol_dist_metabolite * body_mass);
}

void LSTMDrugConversion::setKillingParameters(Params_convFactor& p,
        uint32_t genotype) const
{
    const LSTMDrugPD& pdP = parentType.getPD(genotype), &pdM =
            metaboliteType.getPD(genotype);
    p.nP = pdP.slope();
    p.VP = pdP.max_killing_rate();
    p.KnP = pdP.IC50_pow_slope();
    p.nM = pdM.slope();
    p.VM = pdM.max_killing_rate();
    p.KnM = pdM.IC50_pow_slope();
}

// TODO: in high transmission, is this going to get called more often than updateConcentration?
// When does it make sense to try to optimise (avoid doing decay calcuations here)?
double LSTMDrugConversion::calculateDrugFactor(uint32_t genotype, double body_mass) const {
    if( qtyG == 0.0 && qtyP == 0.0 && qtyM == 0.0 && doses.size() == 0 ){
        return 1.0; // nothing to do
    }
    
    Params_convFactor p;
    setConversionParameters(p, body_mass);
    setKillingParameters(p, genotype);
    
    double time = 0.0;  // time since start of day
    double totalFactor = 1.0;   // survival factor for whole day
    
    typedef pair<double,double> TimeConc;
    foreach( const TimeConc& time_conc, doses ){
        // we iterate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( time < time_conc.first ){
                double duration = time_conc.first - time;
                totalFactor *= calculateFactor(p, duration);
                const double expAbsorb = exp(nka * duration), expPLoss = exp(p.nl * duration);
                p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
                p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
                p.qtyG *= expAbsorb;
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add to quantity of drug in gut:
            p.qtyG += time_conc.second;   // units: mg
        }else /*i.e. tomorrow or later*/{
            // ignore
        }
    }
    if( time < 1.0 ){
        totalFactor *= calculateFactor(p, 1.0 - time);
    }
    
    return totalFactor;
}

void LSTMDrugConversion::updateConcentration( double body_mass ){
    if( qtyG == 0.0 && qtyP == 0.0 && qtyM == 0.0 && doses.size() == 0 ){
        return; // nothing to do
    }
    last_bm = body_mass;
    
    Params_convFactor p;
    setConversionParameters(p, body_mass);
    
    double time = 0.0, duration;
    size_t doses_taken = 0;
    typedef pair<double,double> TimeConc;
    foreach( TimeConc& time_conc, doses ){
        // we iteratate through doses in time order (since doses are sorted)
        if( time_conc.first < 1.0 /*i.e. today*/ ){
            if( (duration = time_conc.first - time) > 0.0 ){
                const double expAbsorb = exp(nka * duration), expPLoss = exp(p.nl * duration);
                p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
                p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
                p.qtyG *= expAbsorb;
                time = time_conc.first;
            }else{ assert( time == time_conc.first ); }
            // add to quantity of drug in gut:
            p.qtyG += time_conc.second;   // units: mg
            doses_taken += 1;
        }else /*i.e. tomorrow or later*/{
            time_conc.first -= 1.0;
        }
    }
    if( time < 1.0 ){
        duration = 1.0 - time;
        const double expAbsorb = exp(nka * duration), expPLoss = exp(p.nl * duration);
        p.qtyM = calculateMetaboliteQuantity(p, expAbsorb, expPLoss, duration);
        p.qtyP = calculateParentQuantity(p, expAbsorb, expPLoss);
        p.qtyG *= expAbsorb;
    }
    //NOTE: would be faster if elements were stored in reverse order — though prescribing would probably be slower
    doses.erase(doses.begin(), doses.begin() + doses_taken);
    
    qtyM = p.qtyM;
    qtyP = p.qtyP;
    qtyG = p.qtyG;
    
    util::streamValidate( qtyM );
    if( qtyP < parentType.getNegligibleConcentration() * body_mass * vol_dist &&
            qtyM < metaboliteType.getNegligibleConcentration() * body_mass * vol_dist_metabolite )
    {
        // once negligible, try to optimise so that we don't have to do
        // anything next time step
        qtyG = qtyP = qtyM = 0.0;
    }
}

void LSTMDrugConversion::checkpoint(ostream& stream){
    qtyG & stream;
    qtyP & stream;
    qtyM & stream;
    nka & stream;
    nkP_sample & stream;
    nconv_sample & stream;
    nkM_sample & stream;
}
void LSTMDrugConversion::checkpoint(istream& stream){
    qtyG & stream;
    qtyP & stream;
    qtyM & stream;
    nka & stream;
    nkP_sample & stream;
    nconv_sample & stream;
    nkM_sample & stream;
}

}
}
