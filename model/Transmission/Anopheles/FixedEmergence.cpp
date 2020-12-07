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

#include "Transmission/Anopheles/FixedEmergence.h"
#include "Transmission/Anopheles/MosqTransmission.h"

#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

#include "rotate.h"

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;


// -----  Initialisation of model, done before human warmup  ------

FixedEmergence::FixedEmergence() :
        initNv0FromSv(numeric_limits<double>::signaling_NaN())
{
    quinquennialS_v.assign (SimTime::fromYearsI(5), 0.0);
    mosqEmergeRate.resize (SimTime::oneYear()); // Only needs to be done here if loading from checkpoint
}


// -----  Initialisation of model which is done after creating initial humans  -----

void FixedEmergence::init2( double tsP_A, double tsP_Amu, double tsP_A1, double tsP_Ah, double tsP_df, double tsP_dff, double EIRtoS_v, MosqTransmission& transmission ){
    // -----  Calculate required S_v based on desired EIR  -----
    
    initNv0FromSv = initNvFromSv * (1.0 - tsP_A - tsP_df);

    // We scale FSCoeffic to give us S_v instead of EIR.
    // Log-values: adding log is same as exponentiating, multiplying and taking
    // the log again.
    FSCoeffic[0] += log( EIRtoS_v );
    vectors::expIDFT (forcedS_v, FSCoeffic, EIRRotateAngle);
    //vectors::expIDFT (forcedS_v, FSCoeffic, FSRotateAngle);
    
    transmission.initState ( tsP_A, tsP_Amu, tsP_A1, tsP_Ah, tsP_df, tsP_dff, initNvFromSv, initOvFromSv, forcedS_v );
    
    // Crude estimate of mosqEmergeRate: (1 - P_A(t) - P_df(t)) / (T * ρ_S) * S_T(t)
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);
    
    // All set up to drive simulation from forcedS_v

    scaleFactor = 1.0;
    shiftAngle = FSRotateAngle;
    scaled = false;
    rotated = false;
}


// -----  Initialisation of model which is done after running the human warmup  -----
bool FixedEmergence::initIterate (MosqTransmission& transmission) {
    // Try to match S_v against its predicted value. Don't try with N_v or O_v
    // because the predictions will change - would be chasing a moving target!
    // EIR comes directly from S_v, so should fit after we're done.

    // Compute avgAnnualS_v from quinquennialS_v for fitting 
    vecDay<double> avgAnnualS_v( SimTime::oneYear(), 0.0 );
    for( SimTime i = SimTime::fromYearsI(4); i < SimTime::fromYearsI(5); i += SimTime::oneDay() ){
        avgAnnualS_v[mod_nn(i, SimTime::oneYear())] =
            quinquennialS_v[i];
    }

    double factor = vectors::sum(forcedS_v) / vectors::sum(avgAnnualS_v);
    
    //cout << "check: " << vectors::sum(forcedS_v) << " " << vectors::sum(avgAnnualS_v) << endl;
    //cout << "Pre-calced Sv, dynamic Sv:\t"<<sumAnnualForcedS_v<<'\t'<<vectors::sum(annualS_v)<<endl;
    if (!(factor > 1e-6 && factor < 1e6)) {
        if( factor > 1e6 && vectors::sum(quinquennialS_v) < 1e-3 ){
            throw util::base_exception("Simulated S_v is approx 0 (i.e.\
 mosquitoes are not infectious, before interventions). Simulator cannot handle this; perhaps\
 increase EIR or change the entomology model.", util::Error::VectorFitting);
        }
        if ( vectors::sum(forcedS_v) == 0.0 ) {
            return false;   // no EIR desired: nothing to do
        }
        cerr << "Input S_v for this vector:\t"<<vectors::sum(forcedS_v)<<endl;
        cerr << "Simulated S_v:\t\t\t"<<vectors::sum(quinquennialS_v)/5.0<<endl;
        throw TRACED_EXCEPTION ("factor out of bounds",util::Error::VectorFitting);
    }

    const double LIMIT = 0.1;

    if(fabs(factor - 1.0) > LIMIT)
    {
        scaled = false;
        double factorDiff = (scaleFactor * factor - scaleFactor) * 1.0;
        scaleFactor += factorDiff;
    }
    else
        scaled = true;

    double rAngle = findAngle(EIRRotateAngle, FSCoeffic, avgAnnualS_v);
    shiftAngle += rAngle;
    rotated = true;

    // cout << "EIRRotateAngle: " << EIRRotateAngle << " rAngle = " << rAngle << ", angle = " << shiftAngle << " scalefactor: " << scaleFactor << " , factor: " << factor << endl;

    // Compute forced_sv from the Fourrier Coeffs
    // shiftAngle rotate the vector to correct the offset between simulated and input EIR
    // shiftAngle is the offset between the 
    vectors::expIDFT(mosqEmergeRate, FSCoeffic, -shiftAngle);
    // Scale the vector according to initNv0FromSv to get the mosqEmergerate
    // scaleFactor scales the vector to correct the ratio between simulated and input EIR
    vectors::scale (mosqEmergeRate, scaleFactor * initNv0FromSv);

    transmission.initIterateScale (factor);//scaleFactor);
    //initNvFromSv *= scaleFactor;     //(not currently used)

    return !(scaled && rotated);
    // return (fabs(factor - 1.0) > LIMIT);// || (rAngle > LIMIT * 2*M_PI / sim::stepsPerYear());
}


double FixedEmergence::update( SimTime d0, double nOvipositing, double S_v ){
    // We use time at end of step (i.e. start + 1) in index:
    SimTime d5Year = mod_nn(d0 + SimTime::oneDay(), SimTime::fromYearsI(5));
    quinquennialS_v[d5Year] = S_v;
    
    // Get emergence at start of step:
    SimTime dYear1 = mod_nn(d0, SimTime::oneYear());
    // Simple model: fixed emergence scaled by larviciding
    return mosqEmergeRate[dYear1] * interventionSurvival();
}

void FixedEmergence::checkpoint (istream& stream){ (*this) & stream; }
void FixedEmergence::checkpoint (ostream& stream){ (*this) & stream; }

}
}
}
