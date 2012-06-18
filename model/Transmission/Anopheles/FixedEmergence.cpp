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

#include "Transmission/Anopheles/FixedEmergence.h"
#include "Transmission/Anopheles/MosqTransmission.h"
#include "Transmission/Anopheles/Nv0DelayFitting.h"

#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;


// -----  Initialisation of model, done before human warmup  ------

FixedEmergence::FixedEmergence() {
    quinquennialS_v.assign (TimeStep::fromYears(5).inDays(), 0.0);
    mosqEmergeRate.resize (TimeStep::DAYS_IN_YEAR); // Only needs to be done here if loading from checkpoint
}


// -----  Initialisation of model which is done after creating initial humans  -----

void FixedEmergence::init2( double tsP_A, double tsP_df, double EIRtoS_v, MosqTransmission& transmission ){
    // -----  Calculate required S_v based on desired EIR  -----
    
    initNv0FromSv = initNvFromSv * (1.0 - tsP_A - tsP_df);

    // We scale FSCoeffic to give us S_v instead of EIR.
    // Log-values: adding log is same as exponentiating, multiplying and taking
    // the log again.
    FSCoeffic[0] += log( EIRtoS_v);
    vectors::expIDFT (forcedS_v, FSCoeffic, FSRotateAngle);
    
    transmission.initState ( tsP_A, tsP_df, initNvFromSv, initOvFromSv, forcedS_v );
    
    // Crude estimate of mosqEmergeRate: (1 - P_A(t) - P_df(t)) / (T * ρ_S) * S_T(t)
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);
    
    // All set up to drive simulation from forcedS_v
}


// -----  Initialisation of model which is done after running the human warmup  -----

bool FixedEmergence::initIterate (MosqTransmission& transmission) {
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
    vector<double> avgAnnualS_v( TimeStep::fromYears(1).inDays(), 0.0 );
    for ( int i = 0; i < TimeStep::fromYears(5).inDays(); ++i ) {
        avgAnnualS_v[i % TimeStep::fromYears(1).inDays()] =
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

    const double LIMIT = 0.1;
    return (fabs(factor - 1.0) > LIMIT) ||
           (rAngle > LIMIT * 2*M_PI / TimeStep::stepsPerYear);
}


double FixedEmergence::get( size_t d, size_t dYear1, double nOvipositing ) {
    //TODO: replace emergence with new formula: rho_p * (num pupae one day from emerging)
    // second two lines don't change
    return mosqEmergeRate[dYear1] * larvicidingIneffectiveness;
}

void FixedEmergence::updateStats( size_t d, double tsP_dif, double S_v ){
    size_t d5Year = d % TimeStep::fromYears(5).inDays();
    quinquennialS_v[d5Year] = S_v;
}

void FixedEmergence::checkpoint (istream& stream){ (*this) & stream; }
void FixedEmergence::checkpoint (ostream& stream){ (*this) & stream; }

}
}
}
