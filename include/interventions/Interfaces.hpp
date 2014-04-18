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
#include "Monitoring/Survey.h"
#include <schema/interventions.h>
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

/** Enumeration of all available components, in the order that these should be
 * deployed in within a single intervention. */
namespace Component { enum Type {
    RECRUIT_ONLY,     // selection for a sub-population (without other effects)
    SCREEN,     // screening, e.g. as part of MSAT
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
struct ComponentId{
    explicit inline ComponentId( size_t id ) : id(id) {}
    explicit inline ComponentId( istream& stream ){ id & stream; }
    inline void operator& (istream& stream) { id & stream; }
    inline void operator& (ostream& stream) const{ id & stream; }
    inline bool operator== (const ComponentId that) const{ return id == that.id; }
    inline bool operator< (const ComponentId that) const{ return id < that.id; }
    size_t id;
};
// special value "whole population cohort" :
static ComponentId ComponentId_pop = ComponentId( boost::integer_traits<size_t>::const_max );

/** A description of one component of a human intervention.
 * 
 * Note that one component can have several "actions", but that deployment and
 * decay of these "actions" is usually related.
 * 
 * This is a base class. */
class HumanInterventionComponent {
public:
    virtual ~HumanInterventionComponent() {}
    
    /** Deploy the component to a pre-selected human.
     * 
     * @param human Individual receiving the intervention
     * @param method Channel of deployment (mass, continuous)
     */
    virtual void deploy( Host::Human& human, Deployment::Method method,
        VaccineLimits vaccLimits ) const =0;
    
    /** Get the component identifier. */
    inline ComponentId id()const{ return m_id; }
    
    /** Get the duration. */
    inline TimeStep duration()const{ return m_duration; }
    
    /** Returns the appropriate descriptor from the Component::Type enum.
     * 
     * This is only used a small number of times during setup, so doesn't need
     * to be fast. */
    virtual Component::Type componentType() const=0;
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const =0;
#endif
    
    // Only for use by InterventionManager:
    inline void setExpireAfter( TimeStep duration ){
        m_duration = duration;
    }
    
protected:
    /** Construct (from a derived class).
     * 
     * @param id Component identifier; used as an identifier for cumulative
     *  deployment as well as to match human-specific components to general
     *  parameters (i.e. objects of the class extending this one).
     * @param ctsMeasure Report measure to use for continuous deployment
     * @param timedMeasure Report measure to use for timed deployment
     */
    explicit HumanInterventionComponent(ComponentId id,
                                        Monitoring::ReportMeasureI ctsMeasure,
                                        Monitoring::ReportMeasureI timedMeasure) :
            m_id(id),
            m_measureCts( ctsMeasure ), m_measureTimed( timedMeasure ) {}
    
    /** Trivial helper function to get deployment measure. */
    inline Monitoring::ReportMeasureI reportMeasure( Deployment::Method method )const{
        return (method == Deployment::TIMED) ? m_measureTimed : m_measureCts;
    }
    
private:
    /** Don't copy (this may be possible but shouldn't be needed). */
    HumanInterventionComponent( const HumanInterventionComponent& );
    
    ComponentId m_id;
    TimeStep m_duration;
    Monitoring::ReportMeasureI m_measureCts, m_measureTimed;
};

/** A description of a human intervention (as a list of components). */
class HumanIntervention {
public:
    /** Create from a list of XML elements: <component id="..."/> **/
    explicit HumanIntervention( const xsd::cxx::tree::sequence<scnXml::Component>& componentList );
    
    /** Deploy all components to a pre-selected human. */
    void deploy( Host::Human& human, Deployment::Method method,
        VaccineLimits vaccLimits ) const;
    
#ifdef WITHOUT_BOINC
    void print_details( std::ostream& out )const;
#endif
    
protected:
    // List of pointers to components. Does not manage memory (InterventionManager::humanComponents does that).
    vector<const HumanInterventionComponent*> components;
};

} }
#endif
