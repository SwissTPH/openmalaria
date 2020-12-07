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

#ifndef Hmod_Anopheles_FixedEmergence
#define Hmod_Anopheles_FixedEmergence

#include "Global.h"
#include "Transmission/Anopheles/EmergenceModel.h"
#include "schema/interventions.h"
#include <vector>
#include <limits>

#include "rotate.h"

namespace OM {
namespace Transmission {
namespace Anopheles {

using namespace std;
using namespace OM::util;

// forward declare to avoid circular dependency:
class MosqTransmission;

/** Part of vector anopheles model, giving emergence of adult mosquitoes from
 * water bodies. This model fits annual (periodic) sequence to produce the
 * desired EIR during warmup, then fixes this level of emergence for the rest
 * of the simulation.
 * 
 * Larviciding intervention directly scales the number of mosquitoes emerging
 * by a number, usually in the range [0,1] (but larger than 1 is also valid).
 */
class FixedEmergence : public EmergenceModel
{
public:
    ///@brief Initialisation and destruction
    //@{
    /// Initialise and allocate memory
    FixedEmergence() : initNv0FromSv(numeric_limits<double>::signaling_NaN())
    {
        quinquennialS_v.assign (SimTime::fromYearsI(5), 0.0);
        mosqEmergeRate.resize (SimTime::oneYear()); // Only needs to be done here if loading from checkpoint
    }

    /** Latter part of AnophelesModel::init2.
     *
     * @param tsP_A P_A for this time step.
     * @param tsP_df P_df for this time step.
     * @param tsP_dff P_dff for this time step.
     * @param EIRtoS_v multiplication factor to convert input EIR into required
     * @param transmission reference to MosqTransmission object
     * S_v. */
    virtual void init2( double tsP_A, double tsP_Amu, double tsP_A1, double tsP_Ah, double tsP_df, double tsP_dff, double EIRtoS_v, MosqTransmission& transmission )
    {
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
    
    /** Work out whether another interation is needed for initialisation and if
     * so, make necessary changes.
     *
     * @returns true if another iteration is needed. */
    virtual bool initIterate (MosqTransmission& transmission)
    {
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
    //@}
    
    virtual double update( SimTime d0, double nOvipositing, double S_v )
    {
        // We use time at end of step (i.e. start + 1) in index:
        SimTime d5Year = mod_nn(d0 + SimTime::oneDay(), SimTime::fromYearsI(5));
        quinquennialS_v[d5Year] = S_v;
        
        // Get emergence at start of step:
        SimTime dYear1 = mod_nn(d0, SimTime::oneYear());
        // Simple model: fixed emergence scaled by larviciding
        return mosqEmergeRate[dYear1] * interventionSurvival();
    }
    
    ///@brief Interventions and reporting
    //@{
    double getResAvailability() const {
        return numeric_limits<double>::quiet_NaN();
    }
    double getResRequirements() const {
        return numeric_limits<double>::quiet_NaN();
    }
    //@}
    
protected:
    virtual void checkpoint (istream& stream){ (*this) & stream; }
    virtual void checkpoint (ostream& stream){ (*this) & stream; }
    
private:
    
    /// Checkpointing
    //Note: below comments about what does and doesn't need checkpointing are ignored here.
    template<class S>
    void operator& (S& stream) {
        mosqEmergeRate & stream;
        quinquennialS_v & stream;
        initNv0FromSv & stream;
    }
    
    // -----  parameters (constant after initialisation)  -----
    
    ///@brief Descriptions of transmission, used primarily during warmup
    //@{
    /** Summary of S_v over the last five years, used by vectorInitIterate to
     * calculate scaling factor.
     *
     * Length is 365 * 5. Checkpoint.
     *
     * Units: inoculations. */
    vecDay<double> quinquennialS_v;
    
    /** Conversion factor from forcedS_v to mosqEmergeRate.
     *
     * Should be checkpointed. */
    double initNv0FromSv;       ///< ditto
    //@}
    
    /** Emergence rate of new mosquitoes, for every day of the year (N_v0).
     * 
     * Has annual periodicity: length is 365. First value (index 0) corresponds
     * to first day of year (1st Jan or something else if rebased). In 5-day
     * time-step model values at indecies 0 through 4 are used to calculate the
     * state at time-step 1.
     * 
     * Units: Animals per day.
     *
     * Should be checkpointed. */
    vecDay<double> mosqEmergeRate;
};

}
}
}
#endif
