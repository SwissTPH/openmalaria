/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2013 Swiss Tropical and Public Health Institute 
 * Copyright (C) 2005-2013 Liverpool School Of Tropical Medicine
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
#ifndef OM_INTERVENTIONS_INTERVENTIONS
#define OM_INTERVENTIONS_INTERVENTIONS

#include "Global.h"
#include "Host/ImportedInfections.h"
#include "schema/interventions.h"
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/concept_check.hpp>
#include <bitset>

namespace OM {
    class Population;
    namespace Host {
        class Human;
    }

namespace interventions {
    using ::boost::ptr_vector;

namespace Deployment {
    enum Method {
        TIMED,   // mass distribution campaign
        CTS     // continuous deployment (EPI, etc.)
    };
}

/** Interface for continuous deployment of an intervention. */
class ContinuousDeployment {
public:
    /// Create, passing deployment age
    explicit ContinuousDeployment( const ::scnXml::ContinuousDeployment& elt );
    virtual ~ContinuousDeployment() {}
    
    /// For sorting
    inline bool operator<( const ContinuousDeployment& that )const{
        return this->deployAge < that.deployAge;
    }
    
    /** Apply filters and potentially deploy.
     * 
     * @returns false iff this deployment (and thus all later ones in the
     *  ordered list) happens in the future. */
    bool filterAndDeploy( Host::Human& human, const Population& population ) const;
    
protected:
    /// Deploy to a selected human.
    virtual void deploy( Host::Human& human, const Population& population ) const =0;
    
    TimeStep begin, end;    // first timeStep active and first timeStep no-longer active
    TimeStep deployAge;
    bool cohortOnly;
    double coverage;
};

/** Interface for timed deployment of an intervention. */
class TimedDeployment {
public:
    /// Create, passing time of deployment
    explicit TimedDeployment(TimeStep deploymentTime);
    virtual ~TimedDeployment() {}
    
    inline bool operator< (const TimedDeployment& that) const{
        return this->time < that.time;
    }
    
    virtual void deploy (OM::Population&) =0;
    
    // Read access required in this file; don't really need protection:
    TimeStep time;
};

/** Enumeration of all effects, in the order that these should be deployed in
 * within a single intervention. */
namespace Effect { enum Type {
    MDA,        // mass drug administration
    PEV,        // pre-erythrocytic vaccine
    BSV,        // blood-stage vaccine
    TBV,        // transmission-blocking vaccine
    IPT,        // intermittent preventative treatment
    ITN,        // insecticide treated net
    IRS,        // indoor residual spraying
    GVI,        // generic vector intervention
    COHORT,     // cohort selection
}; }

/** A description of one effect of a human intervention.
 * 
 * Note that one "effect" can have several "actions", but that deployment and
 * decay of these "actions" is usually related.
 * 
 * This is a base class. */
class HumanInterventionEffect {
public:
    /** Deploy the effect to a pre-selected human.
     * 
     * @param human Individual receiving the intervention
     * @param method Channel of deployment (mass, continuous)
     */
    virtual void deploy( Host::Human& human, Deployment::Method method ) const =0;
    
    /** Get the effect index. */
    inline size_t getIndex()const{ return index; }
    
    /** Returns the appropriate descriptor from the EffectType enum.
     * 
     * This is only used a small number of times during setup, so doesn't need
     * to be fast. */
    virtual Effect::Type effectType() const=0;
    
protected:
    /** Construct (from a derived class).
     * 
     * @param index Effect index; used as an identifier for cumulative
     *  deployment as well as to match human-specific components to general
     *  parameters (i.e. objects of the class extending this one).
     */
    explicit HumanInterventionEffect(size_t index) : index(index) {}
    
private:
    /** Don't copy (this may be possible but shouldn't be needed). */
    HumanInterventionEffect( const HumanInterventionEffect& );
    
    size_t index;
};

/** A description of a human intervention (as a list of effects). */
class HumanIntervention {
public:
    /** Add an effect. */
    inline void addEffect( const HumanInterventionEffect *effect ){ effects.push_back( effect ); }
    
    /** Deploy all effects to a pre-selected human. */
    void deploy( Host::Human& human, Deployment::Method method ) const;
    
    /** Sort effects according to a standard order.
     * 
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's effects (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     * 
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    void sortEffects();
    
private:
    // List of pointers to effects. Does not manage memory (InterventionManager::humanEffects does that).
    vector<const HumanInterventionEffect*> effects;
};

/** Management of interventions deployed on a per-timestep basis. */
class InterventionManager {
public:
    /** Read XML descriptions. */
    InterventionManager (const scnXml::Interventions& intervElt, OM::Population& population);
    
    /// Checkpointing
    template<class S>
    void operator& (S& stream) {
        using namespace OM::util::checkpoint;
        // most members are only set from XML,
        // nextTimed varies but is re-set by loadFromCheckpoint
        importedInfections & stream;
    }

    /** Call after loading a checkpoint, passing the intervention-period time.
     * 
     * Serves to replace health-system and EIR where changeHS/changeEIR
     * interventions have been used. */
    void loadFromCheckpoint( OM::Population& population, OM::util::TimeStep interventionTime );
    
    /// Returns true if any cohort selection "intervention" is active
    inline bool cohortEnabled() const{
        return _cohortEnabled;
    }
    
    /** @brief Deploy interventions
     *
     * Timed interventions are deployed for this timestep.
     * 
     * Continuous interventions are deployed as humans reach the target ages.
     * Unlike with vaccines, missing one schedule doesn't preclude the next. */
    void deploy (OM::Population& population);
    
    /** Get a constant reference to an effect class with a certain index.
     * 
     * @throws util::base_exception if the index is out-of-range */
    inline const HumanInterventionEffect& getEffect( size_t index ){
        if( index >= humanEffects.size() )
            throw util::base_exception( "invalid index" );
        return humanEffects[index];
    }
    
private:
    // All human intervention effects, indexed by a number. This list is used
    // during initialisation and thereafter only for memory management.
    boost::ptr_vector<HumanInterventionEffect> humanEffects;
    // All human interventions. These are stored here for memory management
    // only (so that they are deleted when this class is destroyed).
    boost::ptr_vector<HumanIntervention> humanInterventions;
    // Continuous interventions, sorted by deployment age (weakly increasing)
    ptr_vector<ContinuousDeployment> continuous;
    // List of all timed interventions. Should be sorted (time weakly increasing).
    ptr_vector<TimedDeployment> timed;
    uint32_t nextTimed;
    
    // imported infections are not really interventions, and handled by a separate class
    // (but are grouped here for convenience and due toassociation in schema)
    OM::Host::ImportedInfections importedInfections;
    bool _cohortEnabled;
};

/** Type of container for manager. */
typedef auto_ptr<InterventionManager> auto_ptr_interv_manager;
/** InterventionManager instance. It is the Simulation class's responsibility
 * to set this up before starting simulations. It should exist thereafter. */
extern auto_ptr_interv_manager manager; // defined in Simulation.cpp
} }
#endif
