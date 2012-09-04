/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2012 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Transmission/Anopheles/SimpleMPDEmergence.h"
#include "Transmission/Anopheles/MosqTransmission.h"
#include "Transmission/Anopheles/Nv0DelayFitting.h"

#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "schema/entomology.h"

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;


// -----  Initialisation of model, done before human warmup  ------

SimpleMPDEmergence::SimpleMPDEmergence(const scnXml::SimpleMPD& elt) :
    initNv0FromSv( numeric_limits<double>::quiet_NaN() )
{
    quinquennialS_v.assign (TimeStep::fromYears(5).inDays(), 0.0);
    quinquennialOvipositing.assign (TimeStep::fromYears(5).inDays(), 0.0);
    mosqEmergeRate.assign (TimeStep::DAYS_IN_YEAR, 0.0);
    invLarvalResources.assign (TimeStep::DAYS_IN_YEAR, 0.0);
    
    developmentDuration = elt.getDevelopmentDuration().getValue();
    if (!developmentDuration > 0)
        throw util::xml_scenario_error("entomology.vector.simpleMPD.developmentDuration: "
            "must be positive");
    probPreadultSurvival = elt.getDevelopmentSurvival().getValue();
    if (!(0 <= probPreadultSurvival && probPreadultSurvival <= 1))
        throw util::xml_scenario_error("entomology.vector.simpleMPD.developmentSurvival: "
            "must be a probability (in range [0,1]");
    fEggsLaidByOviposit = elt.getFemaleEggsLaidByOviposit().getValue();
    if (!fEggsLaidByOviposit > 0)
        throw util::xml_scenario_error("entomology.vector.simpleMPD.femaleEggsLaidByOviposit: "
            "must be positive");
    nOvipositingDelayed.assign (developmentDuration, 0.0);
}


// -----  Initialisation of model which is done after creating initial humans  -----

void SimpleMPDEmergence::init2( double tsP_A, double tsP_df, double EIRtoS_v, MosqTransmission& transmission ){
    // -----  Calculate required S_v based on desired EIR  -----
    
    initNv0FromSv = initNvFromSv * (1.0 - tsP_A - tsP_df);

    // We scale FSCoeffic to give us S_v instead of EIR.
    // Log-values: adding log is same as exponentiating, multiplying and taking
    // the log again.
    FSCoeffic[0] += log( EIRtoS_v);
    vectors::expIDFT (forcedS_v, FSCoeffic, FSRotateAngle);
    
    transmission.initState ( tsP_A, tsP_df, initNvFromSv, initOvFromSv, forcedS_v );
    
    // Initialise nOvipositingDelayed
    int y1 = TimeStep::DAYS_IN_YEAR,
        tau = transmission.getMosqRestDuration();
    for (int t=0; t<developmentDuration; ++t){
        nOvipositingDelayed[mod_nn(t+tau, developmentDuration)] =
            tsP_df * initNvFromSv * forcedS_v[t];
    }
    
    // Crude estimate of mosqEmergeRate: (1 - P_A(t) - P_df(t)) / (T * ρ_S) * S_T(t)
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);
    // Used when calculating invLarvalResources (but not a hard constraint):
    assert(tau+developmentDuration <= y1);
    for( int t=0; t<(int)mosqEmergeRate.size(); ++t ){
        double yt = fEggsLaidByOviposit * tsP_df * initNvFromSv *
            forcedS_v[mod_nn(t + y1 - tau - developmentDuration, y1)];
        invLarvalResources[t] = (probPreadultSurvival * yt - mosqEmergeRate[t]) /
            (mosqEmergeRate[t] * yt);
    }
    
    // All set up to drive simulation from forcedS_v
}


// -----  Initialisation of model which is done after running the human warmup  -----

bool SimpleMPDEmergence::initIterate (MosqTransmission& transmission) {
    // Try to match S_v against its predicted value. Don't try with N_v or O_v
    // because the predictions will change - would be chasing a moving target!
    // EIR comes directly from S_v, so should fit after we're done.

    double factor = vectors::sum (forcedS_v)*5 / vectors::sum(quinquennialS_v);
    //cout << "Pre-calced Sv, dynamic Sv:\t"<<sumAnnualForcedS_v<<'\t'<<vectors::sum(annualS_v)<<endl;
    if (!(factor > 1e-6 && factor < 1e6)) {
        if ( vectors::sum(forcedS_v) == 0.0 ) {
            return false;   // no EIR desired: nothing to do
        }
        cerr << "Input S_v for this vector:\t"<<vectors::sum(forcedS_v)<<endl;
        cerr << "Simulated S_v:\t\t\t"<<vectors::sum(quinquennialS_v)/5.0<<endl;
        throw TRACED_EXCEPTION ("factor out of bounds (likely a code error)",util::Error::VectorFitting);
    }

    //cout << "Vector iteration: adjusting with factor "<<factor<<endl;
    // Adjusting mosqEmergeRate is the important bit. The rest should just
    // bring things to a stable state quicker.
    initNv0FromSv *= factor;
    initNvFromSv *= factor;     //(not currently used)
    vectors::scale (mosqEmergeRate, factor);
    transmission.initIterateScale (factor);
    vectors::scale (quinquennialS_v, factor); // scale so we can fit rotation offset

    // average annual period of S_v over 5 years
    vector<double> avgAnnualS_v( TimeStep::DAYS_IN_YEAR, 0.0 );
    for ( int i = 0; i < TimeStep::fromYears(5).inDays(); ++i ) {
        avgAnnualS_v[mod_nn(i, TimeStep::fromYears(1).inDays())] =
            quinquennialS_v[i] / 5.0;
    }

    // Once the amplitude is approximately correct, we try to find a
    // rotation offset.
    double rAngle = Nv0DelayFitting::fit<double> (EIRRotateAngle, FSCoeffic, avgAnnualS_v);
    //cout << "Vector iteration: rotating with angle (in radians): " << rAngle << endl;
    // annualS_v was already rotated by old value of FSRotateAngle, so increment:
    FSRotateAngle -= rAngle;
    vectors::expIDFT (forcedS_v, FSCoeffic, FSRotateAngle);
    // We use the stored initXxFromYy calculated from the ideal population age-structure (at init).
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);
    
    // Finally, update nOvipositingDelayed and invLarvalResources
    vectors::scale (nOvipositingDelayed, factor);
    
    size_t y1 = TimeStep::DAYS_IN_YEAR,
        y2 = 2*TimeStep::DAYS_IN_YEAR,
        y3 = 3*TimeStep::DAYS_IN_YEAR,
        y4 = 4*TimeStep::DAYS_IN_YEAR,
        y5 = 5*TimeStep::DAYS_IN_YEAR;
    assert(mosqEmergeRate.size() == y1);
    
    for( size_t t=0; t<y1; ++t ){
        size_t ttj = t - developmentDuration;
        // b · P_df · avg_N_v(t - θj - τ):
        double yt = fEggsLaidByOviposit * 0.2 * (
            quinquennialOvipositing[ttj + y1] +
            quinquennialOvipositing[ttj + y2] +
            quinquennialOvipositing[ttj + y3] +
            quinquennialOvipositing[ttj + y4] +
            quinquennialOvipositing[mod_nn(ttj + y5, y5)]);
        invLarvalResources[t] = (probPreadultSurvival * yt - mosqEmergeRate[t]) /
            (mosqEmergeRate[t] * yt);
    }
    
    const double LIMIT = 0.1;
    return (fabs(factor - 1.0) > LIMIT) ||
           (rAngle > LIMIT * 2*M_PI / TimeStep::stepsPerYear);
    //NOTE: in theory, mosqEmergeRate and annualEggsLaid aren't needed after convergence.
}


double SimpleMPDEmergence::get( size_t d, size_t dYear1, double nOvipositing ) {
    // Simple Mosquito Population Dynamics model: emergence depends on the
    // adult population, resources available, and larviciding.
    // See: A Simple Periodically-Forced Difference Equation Model for
    // Mosquito Population Dynamics, N. Chitnis, 2012. TODO: publish & link.
    double yt = fEggsLaidByOviposit * nOvipositingDelayed[mod_nn(d, developmentDuration)];
    double emergence = larvicidingIneffectiveness * probPreadultSurvival * yt /
        (1.0 + invLarvalResources[dYear1] * yt);
    nOvipositingDelayed[mod_nn(d, developmentDuration)] = nOvipositing;
    size_t d5Year = mod_nn(d, TimeStep::fromYears(5).inDays());
    quinquennialOvipositing[d5Year] = nOvipositing;
    return emergence;
}

void SimpleMPDEmergence::updateStats( size_t d, double tsP_dif, double S_v ){
    size_t d5Year = mod_nn(d, TimeStep::fromYears(5).inDays());
    quinquennialS_v[d5Year] = S_v;
}

double SimpleMPDEmergence::getResAvailability() const {
    int start = TimeStep::simulation.inDays() - TimeStep::interval;
    double total = 0;
    for (size_t i = 0; i < (size_t)TimeStep::interval; ++i) {
        size_t dYear1 = mod(start + i, TimeStep::DAYS_IN_YEAR);
        total += 1.0 / invLarvalResources[dYear1];
    }
    return total / TimeStep::interval;
}

void SimpleMPDEmergence::checkpoint (istream& stream){ (*this) & stream; }
void SimpleMPDEmergence::checkpoint (ostream& stream){ (*this) & stream; }

}
}
}
