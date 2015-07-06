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

/// Larval resources global data
namespace LR {
    /// Types of resources; the value is the resource type identifier, and the index the key
    vector<size_t> resTypes;
    
    /// For each resource type, a list of species sharing that resource
    vector<vector<size_t> > resUsers;
    
    /** Mean number of female eggs laid when a mosquito oviposites.
     * 
     * Index: species. */
    vector<double> fEggsLaidByOviposit;
    
    /** Duration of development (time from egg laying to emergence) in days.
     * 
     * Index: species. */
    vector<SimTime> developmentDuration;
    
    /** Resources for mosquito larvae (or rather 1 over resources); γ(t) in
     * model description.
     * 
     * Unlike model description, we allow special values 0 for no density
     * dependence and infinity for zero emergence.
     * 
     * The first index t should correspond to the resources available to mosquitoes
     * emerging at time t (i.e. same index in mosqEmergeRate). The second index is the species.
     * 
     * Has annual periodicity: length is 365. First value (index 0) corresponds
     * to first day of year (1st Jan or something else if rebased). In 5-day
     * time-step model values at indecies 0 through 4 are used to calculate the
     * state at time-step 1.
     * 
     * Units: 1 / animals per day.
     *
     * Should be checkpointed. */
    vecDay2D<double> invLarvalResources;
    
    /** Vector for storing values of nOvipositing for the last
     * developmentDuration time steps. Second index 0 should correspond to
     * nOvipositing developmentDuration days before
     * get(0, dYear1, nOvipositing) is called.
     * 
     * First index is species. */
    vector<vecDay<double> > nOvipositingDelayed;
}

void SimpleMPDEmergence::staticCheckpoint(size_t species, ostream& stream){
    LR::invLarvalResources.size2() & stream;
    for( SimTime t = sim::zero(); t < sim::oneYear(); t += sim::oneDay() ){
        LR::invLarvalResources.at(t, species) & stream;
        LR::nOvipositingDelayed[species] & stream;
    }
}
void SimpleMPDEmergence::staticCheckpoint(size_t species, istream& stream){
    size_t dim2;
    dim2 & stream;
    if( LR::invLarvalResources.size2() == 0 ){
        LR::invLarvalResources.assign (sim::oneYear(), dim2, 0.0);
    }
    for( SimTime t = sim::zero(); t < sim::oneYear(); t += sim::oneDay() ){
        LR::invLarvalResources.at(t, species) & stream;
        LR::nOvipositingDelayed[species] & stream;
    }
}


// -----  Initialisation of model, done before human warmup  ------

SimpleMPDEmergence::SimpleMPDEmergence(const scnXml::SimpleMPD& elt, size_t species) :
    species(species),
    initNv0FromSv( numeric_limits<double>::quiet_NaN() )
{
    quinquennialS_v.assign (sim::fromYearsI(5), 0.0);
    quinquennialOvipositing.assign (sim::fromYearsI(5), 0.0);
    mosqEmergeRate.assign (sim::oneYear(), 0.0);
    
    if( LR::developmentDuration.size() <= species ){
        size_t newDim = (species == 0) ? 4 : species * 2;
        LR::developmentDuration.resize( newDim, sim::zero() );
        LR::nOvipositingDelayed.resize( newDim );
        LR::fEggsLaidByOviposit.resize( newDim, numeric_limits<double>::quiet_NaN() );
        LR::invLarvalResources.resize( sim::oneYear(), newDim );
    }
    
    LR::developmentDuration[species] = sim::fromDays(elt.getDevelopmentDuration().getValue());
    if (!(LR::developmentDuration[species] > sim::zero()))
        throw util::xml_scenario_error("entomology.vector.simpleMPD.developmentDuration: "
            "must be positive");
    probPreadultSurvival = elt.getDevelopmentSurvival().getValue();
    if (!(0.0 <= probPreadultSurvival && probPreadultSurvival <= 1.0))
        throw util::xml_scenario_error("entomology.vector.simpleMPD.developmentSurvival: "
            "must be a probability (in range [0,1]");
    
    LR::fEggsLaidByOviposit[species] = elt.getFemaleEggsLaidByOviposit().getValue();
    if (!(LR::fEggsLaidByOviposit[species] > 0.0))
        throw util::xml_scenario_error("entomology.vector.simpleMPD.femaleEggsLaidByOviposit: "
            "must be positive");
    
    LR::nOvipositingDelayed[species].assign( LR::developmentDuration[species], 0.0 );
    
    
    // Resource type identifier. If not present, use species + 1000000 (unique).
    size_t resTypeKey = species + 1000000;
    if( elt.getResourceType().present() ){
        int64_t k = elt.getResourceType().get().getValue();
        if( k < 0 || k >= 1000000 )
            throw util::xml_scenario_error( "entomology.vector.simpleMPD.resourceType must be between 0 and 999999" );
        resTypeKey = k;
    }
    for( resType = 0; ; ++resType ){
        if( resType < LR::resTypes.size() ){
            if( LR::resTypes[resType] == resTypeKey )
                break;      // have our resType
        }else{
            LR::resTypes.push_back( resTypeKey );  // add new resource type
            LR::resUsers.resize( LR::resTypes.size() );
            break;      // resType is LR::resTypes.size()
        }
    }
    LR::resUsers[resType].push_back( species );
}


// -----  Initialisation of model which is done after creating initial humans  -----

void SimpleMPDEmergence::init2( double tsP_A, double tsP_df, double EIRtoS_v, MosqTransmission& transmission ){
    // -----  Calculate required S_v based on desired EIR  -----
    const SimTime devDur = LR::developmentDuration[species];
    
    initNv0FromSv = initNvFromSv * (1.0 - tsP_A - tsP_df);

    // We scale FSCoeffic to give us S_v instead of EIR.
    // Log-values: adding log is same as exponentiating, multiplying and taking
    // the log again.
    FSCoeffic[0] += log( EIRtoS_v);
    vectors::expIDFT (forcedS_v, FSCoeffic, FSRotateAngle);
    
    transmission.initState ( tsP_A, tsP_df, initNvFromSv, initOvFromSv, forcedS_v );
    
    // Initialise nOvipositingDelayed
    SimTime y1 = sim::oneYear();
    SimTime tau = transmission.getMosqRestDuration();
    for( SimTime t = sim::zero(); t < devDur; t += sim::oneDay() ){
        LR::nOvipositingDelayed[species][mod_nn(t+tau, devDur)] =
            tsP_df * initNvFromSv * forcedS_v[t];
    }
    
    // Crude estimate of mosqEmergeRate: (1 - P_A(t) - P_df(t)) / (T * ρ_S) * S_T(t)
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);
    // Used when calculating invLarvalResources (but not a hard constraint):
    assert(tau + devDur <= y1);
    
    for( SimTime t = sim::zero(); t < sim::oneYear(); t += sim::oneDay() ){
        const double y = LR::fEggsLaidByOviposit[species] * tsP_df * initNvFromSv *
            forcedS_v[mod_nn(t + y1 - tau - devDur, y1)];
        const double lr = (probPreadultSurvival * y - mosqEmergeRate[t]) /
            (mosqEmergeRate[t] * y);
        if( lr < 1e-6 ){
            std::cerr << "Low larval resources: " << lr << endl;
            assert( false );
        }
        LR::invLarvalResources.at(t, species) = lr;
    }
    
    // All set up to drive simulation from forcedS_v
}


// -----  Initialisation of model which is done after running the human warmup  -----

bool SimpleMPDEmergence::initIterate (MosqTransmission& transmission) {
    const SimTime devDur = LR::developmentDuration[species];
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
        throw TRACED_EXCEPTION ("factor out of bounds",util::Error::VectorFitting);
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
    vector<double> avgAnnualS_v( sim::oneYear().inDays(), 0.0 );
    for( SimTime i = sim::zero(), end = sim::fromYearsI(5); i < end; i += sim::oneDay() ) {
        avgAnnualS_v[mod_nn(i, sim::oneYear()).inDays()] =
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
    vectors::scale (LR::nOvipositingDelayed[species], factor);
    
    SimTime y1 = sim::oneYear(),
        y2 = sim::fromYearsI(2),
        y3 = sim::fromYearsI(3),
        y4 = sim::fromYearsI(4),
        y5 = sim::fromYearsI(5);
    assert(mosqEmergeRate.size() == y1);
    
    for( SimTime t = sim::zero(); t < y1; t += sim::oneDay() ){
        SimTime ttj = t - devDur;
        // b · P_df · avg_N_v(t - θj - τ):
        const double y = LR::fEggsLaidByOviposit[species] * 0.2 * (
            quinquennialOvipositing[ttj + y1] +
            quinquennialOvipositing[ttj + y2] +
            quinquennialOvipositing[ttj + y3] +
            quinquennialOvipositing[ttj + y4] +
            quinquennialOvipositing[mod_nn(ttj + y5, y5)]);
        const double lr = (probPreadultSurvival * y - mosqEmergeRate[t]) /
            (mosqEmergeRate[t] * y);
        if( lr < 1e-6 ){
            std::cerr << "Low larval resources: " << lr << endl;
            assert( false );
        }
        LR::invLarvalResources.at(t, species) = lr;
    }
    
    const double LIMIT = 0.1;
    return (fabs(factor - 1.0) > LIMIT) ||
           (rAngle > LIMIT * 2*M_PI / sim::stepsPerYear());
    //NOTE: in theory, mosqEmergeRate and annualEggsLaid aren't needed after convergence.
}


double SimpleMPDEmergence::update( SimTime d0, double nOvipositing, double S_v ){
    const SimTime devDur = LR::developmentDuration[species];
    const SimTime d1 = d0 + sim::oneDay();
    const SimTime d5Year = mod_nn(d1, sim::fromYearsI(5));
    quinquennialS_v[d5Year] = S_v;
    
    // Simple Mosquito Population Dynamics model: emergence depends on the
    // adult population, resources available, and larviciding.
    // See: A Simple Periodically-Forced Difference Equation Model for
    // Mosquito Population Dynamics, N. Chitnis, 2012. TODO: publish & link.
    // Modified to allow larval-stage competition.
    
    double y = 0.0, ygamma = 0.0;
    for( vector<size_t>::const_iterator it = LR::resUsers[resType].begin(),
        end = LR::resUsers[resType].end(); it != end; ++it )
    {
        size_t spec = *it;
        double fEggsLaid = LR::fEggsLaidByOviposit[spec] *
            LR::nOvipositingDelayed[spec][mod_nn(d1, devDur)];
        if( spec == species ){
            y = fEggsLaid;
        }
        ygamma += fEggsLaid * LR::invLarvalResources.at(mod_nn(d0, sim::oneYear()), spec);
    }
    
    double emergence = interventionSurvival() * probPreadultSurvival * y /
        (1.0 + ygamma);
    
    if( emergence > 1e6 ){      // for debugging
        std::cerr << "Large emergence!\ny\t" << y << "\nρ\t" << interventionSurvival() * probPreadultSurvival
            << "\ny/γ\t" << ygamma << endl;
        assert(false);
    }
    
    LR::nOvipositingDelayed[species][mod_nn(d1, devDur)] = nOvipositing;
    quinquennialOvipositing[mod_nn(d1, sim::fromYearsI(5))] = nOvipositing;
    return emergence;
}

double SimpleMPDEmergence::getResAvailability() const {
    //TODO: why offset by one time step? This is effectively getting the resources available on the last time step
    //TODO: only have to add one year because of offset
    SimTime start = sim::now() - sim::oneTS() + sim::oneYear();
    double total = 0;
    for( SimTime i = start, end = start + sim::oneTS(); i < end; i += sim::oneDay() ){
        SimTime dYear1 = mod_nn(i, sim::oneYear());
        total += 1.0 / LR::invLarvalResources.at(dYear1, species);
    }
    // Note: these measured resources may be shared; if so the number reported still roughly
    // corresponds to the species.
    return total / sim::oneTS().inDays();
}

void SimpleMPDEmergence::checkpoint (istream& stream){ (*this) & stream; }
void SimpleMPDEmergence::checkpoint (ostream& stream){ (*this) & stream; }

}
}
}
