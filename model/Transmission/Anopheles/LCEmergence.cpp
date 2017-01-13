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

#include "Transmission/Anopheles/LCEmergence.h"
#include "Transmission/Anopheles/MosqTransmission.h"
#include "Transmission/Anopheles/ResourceFitter.h"

#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

#include <cmath>

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;


// -----  Initialisation of model, done before human warmup  ------

LCEmergence::LCEmergence() :
            initialP_A(numeric_limits<double>::quiet_NaN()),
            initialP_df(numeric_limits<double>::quiet_NaN())
{
    quinquennialP_dif.assign (SimTime::fromYearsI(5), 0.0);
#if 0
    mosqEmergeRate.resize (SimTime::oneYear()); // Only needs to be done here if loading from checkpoint
#endif
}

void LCEmergence::initLifeCycle(const scnXml::LifeCycle& lcData){
    lcParams.initLifeCycle( lcData );
    lifeCycle.init( lcParams );
}


// -----  Initialisation of model which is done after creating initial humans  -----

void LCEmergence::init2( double tsP_A, double tsP_df, double EIRtoS_v, MosqTransmission& transmission ){
    initialP_A = tsP_A;
    initialP_df = tsP_df;
    
    // -----  Calculate required S_v based on desired EIR  -----
    
    initNv0FromSv = initNvFromSv * (1.0 - tsP_A - tsP_df);

    // We scale FSCoeffic to give us S_v instead of EIR.
    // Log-values: adding log is same as exponentiating, multiplying and taking
    // the log again.
    FSCoeffic[0] += log( EIRtoS_v);
    vectors::expIDFT (forcedS_v, FSCoeffic, FSRotateAngle);
    
    transmission.initState ( tsP_A, tsP_df, initNvFromSv, initOvFromSv, forcedS_v );
    
    //TODO: VLC merge
    //NOTE: do we not want this first part at all? Or still need some of this?
#if 0
    // Crude estimate of mosqEmergeRate: (1 - P_A(t) - P_df(t)) / (T * ρ_S) * S_T(t)
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);
    // Basic estimate of larvalResources from mosqEmergeRate and human state
    lcParams.fitLarvalResourcesFromEmergence( lcModel, tsP_df, tsP_A, 1, mosqRestDuration, N_v, mosqEmergeRate );
    
    // All set up to drive simulation from forcedS_v
#endif
    
    
    //NOTE: this is only here for debugging. Execution should be the same as
    // that in vectorInitIterate, but much faster.
    cerr << "Warning: LCEmergence::init2: running test mode" << endl;
    const double fixedP_dif[365] = {
        0.0208892, 0.0211708, 0.0211384, 0.0207101, 0.020627, 0.020583,
        0.0204871, 0.0202843, 0.0202156, 0.0200973, 0.0199916, 0.0199837,
        0.0199941, 0.0198548, 0.0197392, 0.0196689, 0.0196281, 0.0196362,
        0.0195606, 0.0196321, 0.0196513, 0.0196813, 0.0197236, 0.019749,
        0.0198011, 0.0199042, 0.0199767, 0.0200928, 0.0203577, 0.020535,
        0.0206487, 0.0207777, 0.020979, 0.0211049, 0.0212751, 0.0213585,
        0.0213085, 0.0215146, 0.0215828, 0.021681, 0.0217439, 0.0218348,
        0.0218485, 0.0218434, 0.0219025, 0.0219072, 0.0218668, 0.0218536,
        0.0218229, 0.0217871, 0.0217856, 0.0217375, 0.0216027, 0.0214194,
        0.0211886, 0.02102, 0.0210537, 0.0210596, 0.020977, 0.0208403,
        0.0209002, 0.0209218, 0.0208807, 0.020998, 0.0209654, 0.0209916,
        0.0210443, 0.0211605, 0.0211625, 0.0211092, 0.0210826, 0.0210332,
        0.0210576, 0.0211502, 0.0212008, 0.0209904, 0.0210191, 0.0210775,
        0.0211045, 0.0209879, 0.0211051, 0.0210568, 0.0209078, 0.0206866,
        0.0204754, 0.0204566, 0.0204005, 0.0204629, 0.0202649, 0.0201598,
        0.0202263, 0.0202866, 0.0202853, 0.0201766, 0.0201573, 0.0201199,
        0.0201266, 0.0202045, 0.0203559, 0.0205059, 0.0205846, 0.020627,
        0.0205889, 0.0207248, 0.0210338, 0.020848, 0.0208437, 0.0208103,
        0.0208386, 0.0210152, 0.0210037, 0.0209833, 0.0209339, 0.0211464,
        0.0211956, 0.021223, 0.0212042, 0.0212204, 0.0212783, 0.0211757,
        0.0211111, 0.0210839, 0.0210912, 0.0210711, 0.0210395, 0.0209084,
        0.0208914, 0.0208499, 0.0208768, 0.0208232, 0.0207477, 0.0206863,
        0.0205907, 0.0204855, 0.0203648, 0.0202348, 0.0201752, 0.0201825,
        0.0202359, 0.0197321, 0.0197388, 0.0197846, 0.019871, 0.0199399,
        0.0199271, 0.0199092, 0.0199943, 0.0199323, 0.019925, 0.0199468,
        0.0198981, 0.0197694, 0.0196401, 0.0196609, 0.0197674, 0.0199079,
        0.0198764, 0.0197805, 0.0196566, 0.019635, 0.0196172, 0.0196139,
        0.0195516, 0.0194965, 0.0194379, 0.0193647, 0.0194341, 0.0195323,
        0.0194786, 0.0194654, 0.0194412, 0.0195673, 0.0196621, 0.0200041,
        0.0201452, 0.0203318, 0.0207462, 0.020807, 0.0208663, 0.0211534,
        0.0214383, 0.0215592, 0.0217233, 0.0218282, 0.0218163, 0.0218591,
        0.0218237, 0.021759, 0.0217767, 0.0217857, 0.0218089, 0.0217847,
        0.0217291, 0.021772, 0.0217833, 0.0217315, 0.0216713, 0.0215465,
        0.021455, 0.0214184, 0.0213807, 0.0213498, 0.0214956, 0.0213265,
        0.0211701, 0.0211106, 0.0212494, 0.0211715, 0.0211351, 0.0211105,
        0.0210394, 0.0210233, 0.0210684, 0.0210343, 0.0210338, 0.0211122,
        0.0211723, 0.0212455, 0.0213314, 0.021379, 0.0213257, 0.0212567,
        0.021173, 0.0210044, 0.020943, 0.0209889, 0.021074, 0.0211136,
        0.0211298, 0.0210848, 0.0208255, 0.0207487, 0.0205689, 0.0203249,
        0.0201748, 0.0200892, 0.0199557, 0.0199712, 0.0200732, 0.0200941,
        0.0200531, 0.0202368, 0.0202247, 0.0202549, 0.0204321, 0.0207082,
        0.0209148, 0.0211646, 0.0214232, 0.0217762, 0.0220687, 0.0222249,
        0.0224237, 0.0223453, 0.0222845, 0.0222731, 0.0223779, 0.0225421,
        0.0226999, 0.0227258, 0.0226618, 0.0225859, 0.0224936, 0.0222001,
        0.0221718, 0.0221996, 0.0221581, 0.0221579, 0.0222042, 0.0221501,
        0.022077, 0.0219974, 0.021962, 0.0218739, 0.0217455, 0.0216513,
        0.0215163, 0.0213792, 0.0213222, 0.0212691, 0.0211987, 0.0211547,
        0.0210673, 0.0210157, 0.0211074, 0.0212388, 0.0213335, 0.0213757,
        0.0215502, 0.0216351, 0.0217362, 0.0217988, 0.0217344, 0.0217221,
        0.0216904, 0.0218072, 0.0219162, 0.0219047, 0.0218265, 0.0217598,
        0.0216885, 0.0216105, 0.0214485, 0.0212917, 0.0211206, 0.0207932,
        0.0207192, 0.0206767, 0.0205668, 0.0204963, 0.0204468, 0.0204232,
        0.0206817, 0.0205926, 0.0206728, 0.0206639, 0.0207806, 0.0208969,
        0.0210323, 0.0212097, 0.0209803, 0.0211011, 0.0212169, 0.0213352,
        0.0215127, 0.0216621, 0.0217956, 0.0218957, 0.0219924, 0.0220033,
        0.021993, 0.0220259, 0.0220868, 0.0222158, 0.0221891, 0.0222084,
        0.0221308, 0.0219296, 0.0217949, 0.0217846, 0.0218326, 0.0219313,
        0.0219141, 0.0218403, 0.0217735, 0.0211396, 0.0210943, 0.021037,
        0.0209607, 0.0209538, 0.0209682, 0.020957, 0.0209123, 0.0209004,
        0.0207769, 0.0206586, 0.0206507, 0.0206338, 0.0206661, 0.0207139,
        0.0206755, 0.0208063, 0.0208376, 0.0209504, 0.0208902 };
    vector<double> fixedP_difVec( fixedP_dif, fixedP_dif + 365 );
    ResourceFitter clm( transmission, lcParams, tsP_A, tsP_df, initNvFromSv, initOvFromSv );
    clm.targetS_vWithP_dif( forcedS_v, fixedP_difVec );
    gsl_vector *x = gsl_vector_alloc( 1 );
    gsl_vector_set_all( x, 1e8 );
    clm.sampler( x );
    //clm.fit();
    exit(16);
}


// -----  Initialisation of model which is done after running the human warmup  -----

bool LCEmergence::initIterate (MosqTransmission& transmission) {
    cerr << "Warning: LCEmergence::initIterate not yet written!" << endl;
    // We now know/can get approximate values for:
    // * human-vector interaction (P_df, P_A) (calculated in init2)
    // * human infectiousness (P_dif) (needs to be sampled over year)
    // * the value of S_v we want to fit to (forcedS_v)
    
    // Find suitable larvalResources using tsP_df and tsP_A
    // TODO: initialise with guessed values for N_v, O_v and S_v
    
    ResourceFitter clm( transmission, lcParams, initialP_A, initialP_df, initNvFromSv, initOvFromSv );
    clm.targetS_vWithP_dif( forcedS_v, quinquennialP_dif );
    clm.fit();
    
    //TODO: free mem? Currently can't because it's still written to.
    //quinquennialP_dif.clear();
    
    // FIXME: run warmup and check resultant EIR
    throw TRACED_EXCEPTION_DEFAULT("TODO");
    return false;
}


double LCEmergence::get( SimTime d0, SimTime dYear1, double nOvipositing ) {
    double emergence = lifeCycle.updateEmergence(lcParams, nOvipositing, d0);
    //TODO
    return emergence * larvicidingIneffectiveness;
}

void LCEmergence::updateStats( SimTime d1, double tsP_dif, double S_v ){
    quinquennialP_dif[mod_nn(d1, SimTime::fromYearsI(5))] = tsP_dif;
}

void LCEmergence::checkpoint (istream& stream){ (*this) & stream; }
void LCEmergence::checkpoint (ostream& stream){ (*this) & stream; }

}
}
}
