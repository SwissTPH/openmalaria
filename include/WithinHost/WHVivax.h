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

#ifndef Hmod_WithinHost_Vivax
#define Hmod_WithinHost_Vivax

#include "Global.h"
#include "WithinHost/WHInterface.h"

#include <list>
#include <memory>

using namespace std;

class UnittestUtil;

namespace scnXml{
    class Primaquine;
}
namespace OM {
namespace WithinHost {
namespace Pathogenesis {
    class PathogenesisModel;
}

/**
 * A brood is the set of hypnozoites resulting from an innoculation, plus a
 * combined blood stage.
 * 
 * In this model, if a hypnozoite releases while a blood stage infection
 * initiated by another hypnozoite from the same brood is active, the newly
 * released hypnozoite does nothing, however, blood stage infections from other
 * broods have no effect.
 */
class VivaxBrood{
public:
    VivaxBrood();
    
    /**
     * Do per-timestep update: remove finished blood stage infections and act
     * on newly releasing hypnozoites.
     * 
     * @param anyNewBloodStage  Set true if any hypnozoite release starts a
     *  blood stage
     * @return true if infection is finished (no more blood or liver stages)
     */
    bool update( bool& anyNewBloodStage );
    
    /** Equivalent to a blood stage existing. We do not model incidence of
     * gametocytes independently, thus this also tests existance of
     * gametocytes. */
    inline bool isPatent() const{
        return bloodStageClearDate > TimeStep::simulation;
    }
    
    /** Fully clear blood stage parasites. */
    void treatmentBS();
    
    /** Fully clear liver stage parasites. */
    void treatmentLS();
    
private:
    // list of times at which the merozoite and hypnozoites release, ordered by
    // time of release, soonest last (i.e. last element is next one to release)
    vector<TimeStep> releaseDates;
    
    // Either TimeStep::never (no blood stage) or a date at which the blood stage will clear.
    TimeStep bloodStageClearDate;
};

/**
 * Implementation of a very basic vivax model.
 * 
 * This is for tropical Vivax (low transmission) and where there is little
 * immunity.
 */
class WHVivax : public WHInterface {
public:
    /// @brief Static methods
    //@{
    /// Initialise static parameters
    static void init(const OM::Parameters& parameters, const scnXml::Scenario& scenario);
    
    /** Set health system parameters (stored in this class for convenience). */
    static void setHSParameters( const scnXml::Primaquine& );
    //@}

    /// @brief Constructors, destructors and checkpointing functions
    //@{
    WHVivax( double comorbidityFactor );
    virtual ~WHVivax();
    //@}
    
    virtual double probTransmissionToMosquito( TimeStep ageTimeSteps, double tbvFactor ) const;
    
    virtual bool summarize(const Host::Human& human);
    
    virtual void importInfection();
    
    virtual void update(int nNewInfs, double ageInYears, double bsvFactor);
    
    virtual bool diagnosticDefault() const;

    virtual Pathogenesis::StatePair determineMorbidity( double ageYears );
    
    virtual void clearImmunity();
    
protected:
    virtual InfectionCount countInfections () const;
    virtual void treatment( TreatmentId treatment );
    virtual bool optionalPqTreatment();
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
private:
    /* Is flagged as never getting PQ: this is a heteogeneity factor. Example:
     * Set to zero if everyone can get PQ, 0.5 if females can't get PQ and
     * males aren't tested (i.e. all can get it) or (1+p)/2 where p is the
     * chance a male being tested and found to be G6PD deficient. */
    bool noPQ;
    
    list<VivaxBrood> infections;
    
    Pathogenesis::State morbidity;

    friend class ::UnittestUtil;
};

}
}
#endif
