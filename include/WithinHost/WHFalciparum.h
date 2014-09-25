/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_WithinHost_Immunity
#define Hmod_WithinHost_Immunity

#include "Global.h"
#include "WithinHost/WHInterface.h"
#include "WithinHost/Treatments.h"

#include <list>
#include <memory>

using namespace std;

class UnittestUtil;

namespace OM {
namespace WithinHost {
namespace Pathogenesis {
    class PathogenesisModel;
}

/**
 * Immunity code and base class for all current P. falciparum models.
 */
class WHFalciparum : public WHInterface {
public:
    /// @brief Static methods
    //@{
    /// Initialise static parameters
    static void init( const OM::Parameters& parameters, const scnXml::Scenario& scenario );
    //@}

    /// @brief Constructors, destructors and checkpointing functions
    //@{
    WHFalciparum( double comorbidityFactor );
    virtual ~WHFalciparum();
    //@}
    
    virtual double probTransmissionToMosquito( double tbvFactor ) const;
    
    virtual bool summarize(const Host::Human& human);
    
    // No PQ treatment for falciparum in current models:
    virtual bool optionalPqTreatment(){ return false; }
    
    virtual inline double getTotalDensity() const{ return totalDensity; }
    
    virtual bool diagnosticDefault() const;
    virtual void treatment( Host::Human& human, TreatmentId treatId );
    virtual void treatSimple(SimTime timeLiver, SimTime timeBlood);
    
    virtual Pathogenesis::StatePair determineMorbidity( double ageYears );

    inline double getCumulative_h() const {
        return m_cumulative_h;
    }
    inline double getCumulative_Y() const {
        return m_cumulative_Y;
    }
    
protected:
    /** Clear infections of the appropriate stages.
     * 
     * @param stage Which stages of the infection are affected
     */
    virtual void clearInfections( Treatments::Stages stage ) =0;
    
    ///@brief Immunity model parameters
    //@{
    /** Updates for the immunity model − assumes m_cumulative_h and m_cumulative_Y
     * have already been incremented.
     *
     * Applies decay of immunity against asexual blood stages, if present. */
    void updateImmuneStatus();

    //!innate ability to control parasite densities
    double _innateImmSurvFact;

    /** Number of infections received since birth. */
    double m_cumulative_h;
    /// Cumulative parasite density since birth (units: days * units of density)
    double m_cumulative_Y;
    /// m_cumulative_Y from previous time step
    double m_cumulative_Y_lag;
    //@}
    
    /// Total asexual blood stage density (sum of density of infections).
    double totalDensity;
    
    /** Maximum parasite density of any infection during the previous interval.
     *
     * With 5-day time steps, this is not just the maximum density of any infection
     * at the end of the time step, but something designed to emulate the maximum
     * of 5 daily samples. */
    double timeStepMaxDensity;
    
    /** Total asexual blood stage density over last 20 days (uses samples from
    * 10, 15 and 20 days ago).
    *
    * m_y_lag[sim::ts0().moduloSteps(y_lag_len)] corresponds to the density
    * from the previous time step (once updateInfection has been called). */
    vector<double> m_y_lag;
    
    /// The PathogenesisModel introduces illness dependant on parasite density
    auto_ptr<Pathogenesis::PathogenesisModel> pathogenesisModel;
    
    /// End of step on which treatment expires = start of first step after expiry
    SimTime treatExpiryLiver, treatExpiryBlood;

    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

    ///@brief Static parameters, set by init()
    //@{
//Standard dev innate immunity for densities
    static double sigma_i;
// contribution of parasite densities to acquired immunity in the presence of fever
    static double immPenalty_22;
    /*
      Remaining immunity against asexual parasites(after time step, each of 2 components y and h)
      This variable decays the effectors m_cumulative_h and m_cumulative_Y in a way that their
      effects on densities (1-Dh and 1-Dy) decay exponentially.
    */
    static double asexImmRemain;
    /*
      Remaining immunity against asexual parasites(after each time step, each of 2 components y and h)
      This variable decays the effectors m_cumulative_h and m_cumulative_Y exponentially.
    */
    static double immEffectorRemain;
    /// Length of m_y_lag array. Wouldn't have to be dynamic if Global::interval was known at compile-time.
    /// set by initHumanParameters
    static int y_lag_len;
    //@}
    
    friend class ::UnittestUtil;
};

}
}
#endif
