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
#ifndef OM_INTERVENTIONS_INTERFACES
#define OM_INTERVENTIONS_INTERFACES

#include "Global.h"
#include <boost/integer_traits.hpp>

namespace scnXml{ class DeploymentBase; }

namespace OM {
    namespace Host {
        class Human;
    }

namespace interventions {

namespace Deployment {
    enum Method {
        TIMED,   // mass distribution campaign
        CTS     // continuous deployment (EPI, etc.)
    };
}

/** Enumeration of all effects, in the order that these should be deployed in
 * within a single intervention. */
namespace Effect { enum Type {
    COHORT,     // cohort selection
    MDA,        // mass drug administration
    MDA_TS1D,   // MDA using the 1-day timestep decision tree and drug action models
    PEV,        // pre-erythrocytic vaccine
    BSV,        // blood-stage vaccine
    TBV,        // transmission-blocking vaccine
    IPT,        // intermittent preventative treatment
    ITN,        // insecticide treated net
    IRS,        // indoor residual spraying
    GVI,        // generic vector intervention
    CLEAR_IMMUNITY,     // reset accumulated immunity to zero
}; }

/** Specifies limits on the number of existing doses when deciding whether to
 * vaccinate a human. */
struct VaccineLimits{
    VaccineLimits() : minPrevDoses( 0 ), maxCumDoses( boost::integer_traits<uint32_t>::const_max ) {}
    void set( const scnXml::DeploymentBase& );
    uint32_t minPrevDoses, maxCumDoses;
};

/** Essentially just an integer, used as a vector index.
 * 
 * This class wraps the integer to guard against unintended conversions to or
 * from an integer (which could be mis-use). */
struct EffectId{
    explicit inline EffectId( size_t id ) : id(id) {}
    explicit inline EffectId( istream& stream ){ id & stream; }
    inline void operator& (istream& stream) { id & stream; }
    inline void operator& (ostream& stream) const{ id & stream; }
    inline bool operator== (const EffectId that) const{ return id == that.id; }
    inline bool operator< (const EffectId that) const{ return id < that.id; }
    size_t id;
};
// special value "whole population cohort" :
static EffectId EffectId_pop = EffectId( boost::integer_traits<size_t>::const_max );

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
    virtual void deploy( Host::Human& human, Deployment::Method method,
        VaccineLimits vaccLimits ) const =0;
    
    /** Get the effect index. */
    inline EffectId id()const{ return m_id; }
    
    /** Returns the appropriate descriptor from the EffectType enum.
     * 
     * This is only used a small number of times during setup, so doesn't need
     * to be fast. */
    virtual Effect::Type effectType() const=0;
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const =0;
#endif
    
protected:
    /** Construct (from a derived class).
     * 
     * @param id Effect index; used as an identifier for cumulative
     *  deployment as well as to match human-specific components to general
     *  parameters (i.e. objects of the class extending this one).
     */
    explicit HumanInterventionEffect(EffectId id) : m_id(id) {}
    
private:
    /** Don't copy (this may be possible but shouldn't be needed). */
    HumanInterventionEffect( const HumanInterventionEffect& );
    
    EffectId m_id;
};

/** A description of a human intervention (as a list of effects). */
class HumanIntervention {
public:
    /** Add an effect. */
    inline void addEffect( const HumanInterventionEffect *effect ){ effects.push_back( effect ); }
    
    /** Deploy all effects to a pre-selected human. */
    void deploy( Host::Human& human, Deployment::Method method,
        VaccineLimits vaccLimits ) const;
    
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
    
#ifdef WITHOUT_BOINC
    void print_details( std::ostream& out )const;
#endif
    
private:
    // List of pointers to effects. Does not manage memory (InterventionManager::humanEffects does that).
    vector<const HumanInterventionEffect*> effects;
};

} }
#endif
