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
#include "WithinHost/Diagnostic.h"
#include "WithinHost/Pathogenesis/State.h"
#include "Parameters.h"

using namespace std;

namespace scnXml{
    class Scenario;
    class Model;
    class TreatmentOption;
}
class UnittestUtil;

namespace OM {
namespace Host {
    class Human;
}
namespace WithinHost {

/**
 * Type used to select a treatment option.
 * 
 * Pass by value; it just hides an integer.
 * 
 * NOTE: this struct and the Treatments class are no longer strictly necessary.
 * Ideally, one would remove this, firstly by replacing usages of
 * om:TreatmentOption in healthSystem.xsd with om:DecisionTree (requires an
 * XML updator algorithm), then optionally removing ImmediateOutcomes.
 */
//TODO: can we construct a CMDecisionTree when loading the relevant bits from the XML instead of changing the XSD?
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
     * This step is independent of parasite genetics.
     * 
     * @param tbvFactor Probability that transmission is not blocked by a
     *  "transmission blocking vaccine".
     * @param sumX Optional out parameter for usage with
     *  pTransGenotype(). Pointer may be zero if not required, and is
     *  expected to be zero when not using WHFalciparum.
     * @returns the probability of this human infecting a feeding mosquito.
     * 
     * Calculates the value during the call, which is expensive (cache externally
     * if the value is needed multiple times). */
    virtual double probTransmissionToMosquito( double tbvFactor,
                                               double *sumX )const =0;
    /** Calculates a probability of transmitting an infection of a given
     * genotype to a mosquito, given the two outputs of
     * probTransmissionToMosquito(). Only available for WHFalciparum and
     * should only be called when num genotypes > 1. */
    inline double probTransGenotype( double pTrans, double sumX, size_t genotype ){
        if( pTrans <= 0.0 ) return 0.0;
        else return pTransGenotype( pTrans, sumX, genotype );
    }
    
    /// @returns true if host has patent parasites
    virtual bool summarize(const Host::Human& human) =0;

    /// Create a new infection within this human
    virtual void importInfection() =0;

    /**
     * Carry out the effects of some treatment option.
     * 
     * This may be used by intervention deployment, thus should use time 
     * sim::nowOrTs1(). */
    //TODO: remove this (use treatSimple instead)
    virtual void treatment( Host::Human& human, TreatmentId treatId ) =0;
    
    /** Conditionally gives Primaquine as a treatment.
     * 
     * Returns true iff PQ is administered. Administered implies either fully
     * effective or no effect, depending on another probability. Not
     * administered implies no effect. */
    virtual bool optionalPqTreatment() =0;
    
    /** Treat a patient via the simple treatment model. */
    virtual void treatSimple(SimTime timeLiver, SimTime timeBlood) =0;
    
    /** Give a patient a course of drugs, via the Pk/Pd model
     * 
     * Note: doses sizes are modified according to age via the dosage
     * table given at the time this function is called.
     *
     * @param schedule Index of a treatment schedule
     * @param dosages Index of a dosage table
     * @param age Age of human in years
     */
    virtual void treatPkPd(size_t schedule, size_t dosages, double age) =0;

    /** Add new infections and update the parasite densities of existing
     * infections. Also update immune status.
     *
     * @param nNewInfs Number of inoculations this time-step
     * @param genotype_weights For use in selecting infection genotypes. See
     *  documentation of Genotypes::sampleGenotype().
     * @param ageInYears Age of human
     * @param bsvFactor Parasite survival factor for blood-stage vaccines
     * @param drugMon Only required for a drug monitoring HACK and could be
     *  removed
     */
    virtual void update(int nNewInfs, vector<double>& genotype_weights,
            double ageInYears, double bsvFactor, ofstream& drugMon) =0;

    /** TODO: this should not need to be exposed. It is currently used by a
     * severe outcome (pDeath) model inside the EventScheduler "case
     * management" model, and case management diagnostics. */
    virtual double getTotalDensity() const =0;
    
    /** Simulate use of a diagnostic test.
     *
     * Does not report for costing purposes.
     * 
     * Is used both during time step updates and during monitoring.
     * 
     * @returns true when the diagnostic is positive
     */
    virtual bool diagnosticResult( const Diagnostic& diagnostic ) const =0;
    
    /** Use the pathogenesis model to determine, based on infection status
     * and random draw, this person't morbidity.
     * 
     * This function is called after update() every time step.
     * 
     * @param ageYears Age of human host in years
     */
    virtual Pathogenesis::StatePair determineMorbidity( double ageYears ) =0;

    /// Special intervention: clears all immunity
    virtual void clearImmunity() =0;
    
    // TODO(monitoring): these shouldn't have to be exposed (perhaps use summarize to report the data):
    virtual double getCumulative_h() const =0;
    virtual double getCumulative_Y() const =0;

    /** The maximum number of infections a human can have. The only real reason
     * for this limit is to prevent incase bad input from causing the number of
     * infections to baloon stupidly.
     *
     * Exact constraint is: _MOI <= MAX_INFECTIONS. */
    static const int MAX_INFECTIONS = 21;

protected:
    // See probTransGenotype; this function should only be called when pTrans > 0
    virtual double pTransGenotype( double pTrans, double sumX,
                                   size_t genotype ) =0;
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);

    /// Multiplicity of infection
    int numInfs;

    friend class ::UnittestUtil;
};

}
}
#endif
