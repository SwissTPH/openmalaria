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
#include "Parameters.h"

using namespace std;

namespace scnXml{
    class Scenario;
    class TreatmentOption;
}
class UnittestUtil;

namespace OM {
namespace WithinHost {

/**
 * Type used to select a treatment option.
 * 
 * Pass by value; it just hides an integer.
 */
struct TreatmentId{
    inline bool operator==( const TreatmentId that ){ return id == that.id; }
    inline bool operator!=( const TreatmentId that ){ return id != that.id; }
    
    /// Default constructor: construct to an initial value. Don't pass this value to WHInterface::treatment()!
    TreatmentId() : id(numeric_limits<uint32_t>::max()) {}
private:
    explicit TreatmentId(uint32_t id) : id(id) {}
    uint32_t id;
    friend class Treatments;
};

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
    static void init(const OM::Parameters& parameters, const scnXml::Scenario& scenario);
    
    /** Configure a new treatment option, and return the code used to select
     * that option later. */
    static TreatmentId addTreatment( const scnXml::TreatmentOption& desc );

    /// Create an instance using the appropriate model
    static WHInterface* createWithinHostModel( double comorbidityFactor );
    //@}

    /// @brief Constructors, destructors and checkpointing functions
    //@{
    WHInterface();
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
    virtual double probTransmissionToMosquito( TimeStep ageTimeSteps, double tbvFactor ) const =0;
    
    /// @returns true if host has patent parasites
    virtual bool summarize(const Host::Human& human) =0;

    /// Create a new infection within this human
    virtual void importInfection() =0;

    /**
     * Carry out the effects of some treatment option
     */
    virtual void treatment( Host::Human& human, TreatmentId treatId ) =0;
    
    /** Conditionally gives Primaquine as a treatment.
     * 
     * Returns true iff PQ is administered. Administered implies either fully
     * effective or no effective, depending on another probability. Not
     * administered implies no effect. */
    virtual bool optionalPqTreatment() =0;
    
    /** Give a patient a course of drugs, via the Pk/Pd model
     * 
     * Note: doses sizes are modified according to age via the dosage
     * table given at the time this function is called.
     *
     * @param schedule Index of a treatment schedule
     * @param dosages Index of a dosage table
     */
    virtual void treatPkPd(size_t schedule, size_t dosages) =0;

    /** Add new infections and update the parasite densities of existing
     * infections. Also update immune status.
     *
     * @param nNewInfs Number of inoculations this time-step
     * @param ageInYears Age of human
     * @param bsvFactor Parasite survival factor for blood-stage vaccines */
    virtual void update(int nNewInfs, double ageInYears, double bsvFactor) =0;

    /** TODO: this should not need to be exposed
     * 
     * It is used by: MDA diagnostics, EventScheduler diagnostics, and a severe
     * outcome (pDeath) model inside the EventScheduler "case management"
     * model. */
    virtual double getTotalDensity() const =0;
    
    /** Simulate use of a diagnostic test, using the general detection limit.
     * Does not report for costing purposes.
     * 
     * @returns true when the diagnostic is positive
     */
    virtual bool diagnosticDefault() const =0;
    
    /** Use the pathogenesis model to determine, based on infection status
     * and random draw, this person't morbidity.
     * 
     * This function is called after update() every timestep.
     * 
     * @param ageYears Age of human host in years
     */
    virtual Pathogenesis::StatePair determineMorbidity( double ageYears ) =0;

    /// Special intervention: clears all immunity
    virtual void clearImmunity() =0;
    
    // TODO(monitoring): these shouldn't have to be exposed (perhaps use summarize to report the data):
    virtual double getCumulativeh() const =0;
    virtual double getCumulativeY() const =0;

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

    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

    /// Multiplicity of infection
    int numInfs;

    friend class ::UnittestUtil;
};

}
}
#endif
