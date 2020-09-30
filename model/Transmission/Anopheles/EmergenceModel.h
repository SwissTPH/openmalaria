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

#ifndef Hmod_Anopheles_EmergenceModel
#define Hmod_Anopheles_EmergenceModel

#include "Global.h"
#include "schema/interventions.h"
#include "util/SimpleDecayingValue.h"
#include "util/checkpoint_containers.h"
#include "util/vecDay.h"

namespace OM {
namespace Transmission {
namespace Anopheles {

using namespace std;
using util::vecDay;
using util::LocalRng;

// forward declare to avoid circular dependency:
class MosqTransmission;

/** Part of vector anopheles model, giving emergence of adult mosquitoes from
 * water bodies.
 * 
 * This is an abstract class (i.e. interface). The following implementations
 * exist: FixedEmergence.
 */
class EmergenceModel
{
public:
    ///@brief Initialisation and destruction
    //@{
    EmergenceModel();
    virtual ~EmergenceModel() {}
    
    /** Called to initialise life-cycle parameters from XML data.
     * 
     * Default implementation throws an assertion: only call if the life cycle
     * model is used. */
    virtual void initLifeCycle(const scnXml::LifeCycle& lcData){
        assert (false);
    }
    
    /** Called by initialise function to init variables directly related to EIR
     * 
     * @param anoph Data from XML
     * @param initialisationEIR In/out parameter: TransmissionModel::initialisationEIR
     * @param EIPDuration parameter from MosqTransmission (used for an estimation)
     */
    void initEIR(
        const scnXml::AnophelesParams& anoph,
        vector<double>& initialisationEIR,
        SimTime EIPDuration );
    
    /** Set up the non-host-specific interventions. */
    void initVectorInterv( const scnXml::VectorSpeciesIntervention& elt, size_t instance );
    
    /** Scale the internal EIR representation by factor; used as part of
     * initialisation. */
    void scaleEIR( double factor );
    
    /** Latter part of AnophelesModel::init2.
     *
     * @param tsP_A P_A for initialisation; assumed constant when there are no interventions
     * @param tsP_df P_df for initialisation; assumed constant when there are no interventions
     * @param tsP_dff P_dff for initialisation; assumed constant when there are no interventions
     * @param EIRtoS_v multiplication factor to convert input EIR into required
     * @param transmission reference to MosqTransmission object
     * S_v. */
    virtual void init2( double tsP_A, double tsP_df, double tsP_dff, double EIRtoS_v, MosqTransmission& transmission ) =0;
    
    /** Work out whether another interation is needed for initialisation and if
     * so, make necessary changes.
     *
     * @returns true if another iteration is needed. */
    //virtual bool initIterate (MosqTransmission& transmission, int populationSize, const vector<double> &laggedkappa) =0;
    //@}
    
    /// Update per time-step (for larviciding intervention). Call before
    /// getting emergence each time-step.
    void update ();
    
    vecDay<double> &getEIR() { return speciesEIR; };

    vecDay<double> &getInputEIR() { return speciesEIR; };
    
    /** Model updates.
     * 
     * Returns the emergence for today, taking interventions like larviciding
     * into account, and updates some statistics (needed during
     * initialisation).
     * 
     * @param d0 Time of the start of the day-long update period
     * @param nOvipositing The number of adults which successfully
     * oviposited this/last time-step.
     * @param S_v Value of S_v for this day
     * @returns The number of adults emerging between the last simulated time
     * point and the one being calculated.
     */
    virtual double update( SimTime d0, double nOvipositing, double S_v ) =0;
    
    virtual vecDay<double> &getEmergenceRate() = 0;

    ///@brief Interventions and reporting
    //@{
    /// Start an intervention affecting the vector population.
    void deployVectorPopInterv (LocalRng& rng, size_t instance);
    
    virtual double getResAvailability() const =0;
    virtual double getResRequirements() const =0;
    //@}

    const double getinitNvFromSv() const { return initNvFromSv; }
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        EIRRotateAngle & stream;
        FSRotateAngle & stream;
        FSCoeffic & stream;
        forcedS_v & stream;
        initNvFromSv & stream;
        initOvFromSv & stream;
        emergenceReduction & stream;
        emergenceSurvival & stream;
        checkpoint (stream);
    }
    
protected:
    /** Return the proportion of emerging larvae to survive intervention
     * effects. Should be between 0 and 1. */
    inline double interventionSurvival() const{ return emergenceSurvival; }
    
    virtual void checkpoint (istream& stream) =0;
    virtual void checkpoint (ostream& stream) =0;
    
    
    // -----  parameters (constant after initialisation)  -----
    
    ///@brief Descriptions of transmission, used primarily during warmup
    //@{
    /// Angle (in radians) to rotate series generated by FSCoeffic by, for EIR.
    double EIRRotateAngle;

    /// Rotation angle (in radians) for emergence rate. Both offset for EIR given in XML file and
    /// offset needed to fit target EIR (delayed from emergence rate). Checkpoint.
    double FSRotateAngle;

    /** Fourier coefficients for EIR / forcedS_v series, input from XML file.
     *
     * Initially used to calculate initialisation EIR, then scaled to calc. S_v.
     *
     * When calcFourierEIR is used to produce an EIR from this over 365
     * (365) elements, the resulting EIR has units of
     * infectious bites per adult per day.
     *
     * fcEir must have odd length and is ordered: [a0, a1, b1, ..., an, bn].
     * FSCoeffic[0] needs checkpointing, the rest doesn't. */
    vector<double> FSCoeffic;

    /** S_v used to force an EIR during vector init.
     * 
     * Has annual periodicity: length is 365. First value (index 0) corresponds
     * to first day of year (1st Jan or something else if rebased). In 5-day
     * time-step model values at indecies 0 through 4 are used to calculate the
     * state at time-step 1.
     *
     * Should be checkpointed. */
    vecDay<double> forcedS_v;
    
    /** Conversion factor from forcedS_v to (initial values of) N_v (1 / ρ_S).
     * Should be checkpointed. */
    double initNvFromSv;
    
    /** Conversion factor from forcedS_v to (initial values of) O_v (ρ_O / ρ_S).
     * Should be checkpointed. */
    double initOvFromSv;
    //@}
    
    /** @brief Intervention parameters
     *
     * Checkpointed. */
    //@{
    /// Description of intervention killing effects on emerging pupae
    vector<util::SimpleDecayingValue> emergenceReduction;
    /// Cache parameter updated by update()
    double emergenceSurvival;   // survival with regards to intervention effects
    //@}

    vecDay<double> speciesEIR;

    double scaleFactor, shiftAngle;
    bool rotated, scaled;
};

}
}
}
#endif
