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

#ifndef Hmod_WithinHost_Interface
#define Hmod_WithinHost_Interface

#include "Global.h"
#include "Monitoring/Survey.h"
#include "WithinHost/Pathogenesis/State.h"

#include <list>

using namespace std;

class UnittestUtil;

namespace OM {
namespace WithinHost {

/**
 * Interface to the within-host models. These models encapsulate the infections
 * and related immunity factors of a single human, starting with infection
 * (i.e. assuming successful innoculation), including some drug action code,
 * and outputting parasite densities.
 */
class WHInterface {
public:
    /// @brief Static methods
    //@{
    /// Initialise static parameters
    static void init();

    /// Create an instance using the appropriate model
    static WHInterface* createWithinHostModel ();
    //@}

    /// @brief Constructors, destructors and checkpointing functions
    //@{
    WHInterface();
    /** Second step of initialisation (could be combined with constructor, but
     * for the moment separate to avoid changing the order of random number
     * samples). */
    virtual void setComorbidityFactor( double factor ) =0;
    virtual ~WHInterface();

    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        checkpoint (stream);
    }
    //@}

    /** Return the infectiousness of this human to biting mosquitoes.
     * 
     * @param ageTimeSteps Age of the human
     * 
     * Calculates the value during the call, which is expensive (cache externally
     * if the value is needed multiple times). */
    //TODO: per genotype? (for LSTM's spread of resistance modelling)
    virtual double probTransmissionToMosquito( TimeStep ageTimeSteps, double tbvEfficacy ) const =0;
    
    /// @returns true if host has patent parasites
    virtual bool summarize(Monitoring::Survey& survey, Monitoring::AgeGroup ageGroup) =0;

    /// Create a new infection within this human
    virtual void importInfection() =0;
    /** Conditionally clears all infections. Not used with the PK/PD model.
     *
     * If IPT isn't present, this literally removes all infections, both blood
     * stage and liver stage (TODO: should it remove liver stage?).
     * 
     * When using the IPT model, this conditionally either does the above or
     * does nothing.
     */
    virtual void clearInfections (bool isSevere);

    /** Medicate drugs (wraps drug's medicate).
     *
     * @param drugAbbrev	abbrevation of drug name (e.g. CQ, MF)
     * @param qty		Quantity of drug to administer in mg
     * @param time		Time relative to beginning of timestep to medicate at, in days (less than 1 day)
     * @param duration Duration in days. 0 or NaN indicate oral treatment.
     * @param bodyMass	Weight of human in kg
     */
    virtual void medicate(string drugAbbrev, double qty, double time, double duration, double bodyMass);

    /** Add new infections and update the parasite densities of existing
     * infections. Also update immune status.
     *
     * @param nNewInfs Number of inoculations this time-step
     * @param ageInYears Age of human
     * @param BSVEfficacy Efficacy of blood-stage vaccine */
    virtual void update(int nNewInfs, double ageInYears, double BSVEfficacy) =0;

    // TODO: these should not be exposed outsite the withinhost models,
    // but must be until the 1-day time step case management models are updated
    virtual double getTotalDensity() const;
    
    /** Simulate use of a diagnostic test, using the general detection limit.
     * Does not report for costing purposes.
     * 
     * @returns true when the diagnostic is positive
     */
    virtual bool diagnosticDefault() const =0;
    
    /** Simulate use of a diagnostic test, using the parameters defined for use
     * with MSAT. Does not report for costing purposes.
     * 
     * @returns true when the diagnostic is positive
     * 
     * TODO: this should be generalised and used for more than just MDA.
     */
    virtual bool diagnosticMDA() const =0;
    
    /** Use the pathogenesis model to determine, based on infection status
     * and random draw, this person't morbidity.
     * 
     * @param ageYears Age of human host in years
     */
    virtual Pathogenesis::State determineMorbidity( double ageYears ) =0;

    ///@brief Only do anything when IPT is present:
    //@{
    /// Continuous deployment for IPT
    virtual void continuousIPT (Monitoring::AgeGroup ageGroup, bool inCohort);
    /// Timed deployment for IPT
    virtual void timedIPT (Monitoring::AgeGroup ageGroup, bool inCohort);
    /// Last IPTi dose recent enough to give protection?
    virtual bool hasIPTiProtection (TimeStep maxInterventionAge) const;
    //@}
    
    /// Special intervention: clears all immunity
    virtual void immuneSuppression() =0;

    // TODO: these shouldn't have to be exposed (perhaps use summarize to report the data):
    virtual double getCumulativeh() const;
    virtual double getCumulativeY() const;

    /** The maximum number of infections a human can have. The only real reason
     * for this limit is to prevent incase bad input from causing the number of
     * infections to baloon stupidly.
     *
     * Exact constraint is: _MOI <= MAX_INFECTIONS. */
    static const int MAX_INFECTIONS = 21;

protected:

    struct InfectionCount{
        InfectionCount(): total(0), patent(0) {}        // initialise to 0
        int total;      // includes blood and liver stages
        int patent;     // number of detectible blood-stage infections
    };
    /** For summarizing:
     * 
     * @returns Number of infections, patent and total
     */
    virtual InfectionCount countInfections () const =0;

    /** Literally just removes all infections in an individual.
     *
     * Normally clearInfections() would be called instead, which, when IPT is not
     * active, just calls this function (although this needs to be changed for
     * PK_PD integration). */
    virtual void clearAllInfections() =0;

    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

    /// Multiplicity of infection
    int numInfs;

    friend class ::UnittestUtil;
};

}
}
#endif
