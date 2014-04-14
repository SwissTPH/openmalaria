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
#include <boost/ptr_container/ptr_list.hpp>

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

class WHVivax;
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
    /** Create.
     * 
     * @param host      The host creating this. Not really needed, except to
     * prevent VivaxBrood being default-constructed by its container.
     */
    VivaxBrood( WHVivax *host );
    ~VivaxBrood();
    /** Save a checkpoint. */
    void checkpoint( ostream& stream );
    /** Create from checkpoint. */
    VivaxBrood( istream& stream );
    
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
    VivaxBrood() {}     // not default constructible
    VivaxBrood( const VivaxBrood& ) {}  // not copy constructible
    
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
    static void init();
    
    /** Set health system parameters (stored in this class for convenience). */
    static void setHSParameters( const scnXml::Primaquine& );
    //@}

    /// @brief Constructors, destructors and checkpointing functions
    //@{
    WHVivax();
    virtual void setComorbidityFactor( double factor );
    virtual ~WHVivax();
    //@}
    
    virtual double probTransmissionToMosquito( TimeStep ageTimeSteps, double tbvEfficacy ) const;
    
    virtual bool summarize(Monitoring::Survey& survey, Monitoring::AgeGroup ageGroup);
    
    virtual void importInfection();
    
    virtual void update(int nNewInfs, double ageInYears, double BSVEfficacy);
    
    virtual bool diagnosticDefault() const;
    virtual bool diagnosticMDA() const;

    virtual Pathogenesis::State determineMorbidity( double ageYears );
    
    virtual void immuneSuppression();
    
protected:
    virtual InfectionCount countInfections () const;
    virtual void effectiveTreatment();
    virtual bool optionalPqTreatment();
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
private:
    WHVivax( const WHVivax& ) {}        // not copy constructible
    
    /* Is flagged as never getting PQ: this is a heteogeneity factor. Example:
     * Set to zero if everyone can get PQ, 0.5 if females can't get PQ and
     * males aren't tested (i.e. all can get it) or (1+p)/2 where p is the
     * chance a male being tested and found to be G6PD deficient. */
    bool noPQ;
    
    //TODO: shouldn't have to store by pointer (at least, if using C++11)
    boost::ptr_list<VivaxBrood> infections;
    
    Pathogenesis::State morbidity;

    friend class ::UnittestUtil;
};

}
}
#endif
