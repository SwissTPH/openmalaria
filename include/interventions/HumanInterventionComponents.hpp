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

// Note: this file is included by exactly one source (Interventions.cpp)
// and contains definitions as well as declarations.

// The includes here are more for documentation than required.
#include "interventions/Interfaces.hpp"
#include "interventions/InterventionManager.hpp"
#include "Monitoring/Survey.h"
#include "util/ModelOptions.h"
#include "Clinical/ESCaseManagement.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Clinical/CaseManagementCommon.h"
#include "Host/Human.h"
#include "Transmission/TransmissionModel.h"
#include "WithinHost/Diagnostic.h"
#include "util/random.h"
#include <schema/healthSystem.h>
#include <schema/interventions.h>

namespace OM { namespace interventions {
    using Host::Human;
    using namespace Monitoring;

// ———  HumanIntervention  ———

bool componentCmp(const HumanInterventionComponent *a, const HumanInterventionComponent *b){
    return a->componentType() < b->componentType();
}
HumanIntervention::HumanIntervention( const xsd::cxx::tree::sequence<scnXml::Component>& componentList ){
    components.reserve( componentList.size() );
    for( xsd::cxx::tree::sequence<scnXml::Component>::const_iterator it = componentList.begin(),
        end = componentList.end(); it != end; ++it )
    {
        ComponentId id = InterventionManager::getComponentId( it->getId() );
        const HumanInterventionComponent* component = &InterventionManager::getComponent( id );
        components.push_back( component );
    }
    
    /** Sort components according to a standard order.
     * 
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's components (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     * 
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    std::stable_sort( components.begin(), components.end(), componentCmp );
}
HumanIntervention::HumanIntervention( const xsd::cxx::tree::sequence<scnXml::DTDeploy>& componentList ){
    components.reserve( componentList.size() );
    for( xsd::cxx::tree::sequence<scnXml::DTDeploy>::const_iterator it = componentList.begin(),
        end = componentList.end(); it != end; ++it )
    {
        ComponentId id = InterventionManager::getComponentId( it->getComponent() );
        const HumanInterventionComponent* component = &InterventionManager::getComponent( id );
        components.push_back( component );
    }
    
    /** Sort components according to a standard order.
     * 
     * The point of this is to make results repeatable even when users change
     * the ordering of a list of intervention's components (since getting
     * repeatable results out of OpenMalaria is often a headache anyway, we
     * might as well at least remove this hurdle).
     * 
     * Note that when multiple interventions are deployed simultaneously, the
     * order of their deployments is still dependent on the order in the XML
     * file. */
    std::stable_sort( components.begin(), components.end(), componentCmp );
}

void HumanIntervention::deploy( Human& human, Deployment::Method method,
    VaccineLimits vaccLimits ) const
{
    for( vector<const HumanInterventionComponent*>::const_iterator it = components.begin();
            it != components.end(); ++it )
    {
        const interventions::HumanInterventionComponent& component = **it;
        // we must report first, since it can change cohort and sub-population
        // which may affect what deployment does (at least in the case of reporting deployments)
        human.reportDeployment( component.id(), sim::fromTS(component.duration()) );
        component.deploy( human, method, vaccLimits );
    }
}

#ifdef WITHOUT_BOINC
void HumanIntervention::print_details( std::ostream& out )const{
    out << "human:";
    for( vector<const HumanInterventionComponent*>::const_iterator it =
        components.begin(); it != components.end(); ++it ){
        out << '\t' << (*it)->id().id;
    }
}
#endif

// ———  Utilities  ———

TriggeredDeployments::TriggeredDeployments(const scnXml::TriggeredDeployments& elt){
    lists.reserve( elt.getDeploy().size() );
    for( scnXml::TriggeredDeployments::DeployConstIterator it = elt.getDeploy().begin(),
        end = elt.getDeploy().end(); it != end; ++it )
    {
        lists.push_back( SubList( *it ) );
    }
}
void TriggeredDeployments::deploy(Human& human,
        Deployment::Method method, VaccineLimits vaccLimits) const
{
    for( vector<SubList>::const_iterator it = lists.begin(); it != lists.end(); ++it ){
        it->deploy( human, method, vaccLimits );
    }
}

TriggeredDeployments::SubList::SubList( const scnXml::TriggeredDeployments::DeployType& elt ) :
        HumanIntervention( elt.getComponent() ),
        minAge( TimeStep::fromYears( elt.getMinAge() ) ),
        maxAge( TimeStep::future ),
        coverage( elt.getP() )
{
    if( elt.getMaxAge().present() )
        maxAge = TimeStep::fromYears( elt.getMaxAge().get() );
    
    if( minAge < TimeStep(0) || maxAge < minAge ){
        throw util::xml_scenario_error("triggered intervention must have 0 <= minAge <= maxAge");
    }
    
    if( coverage < 0.0 || coverage > 1.0 ){
        throw util::xml_scenario_error( "triggered intervention must have 0 <= coverage <= 1" );
    }
    
    // Zero coverage: optimise
    if( coverage <= 0.0 || minAge >= maxAge ) components.clear();
}
void TriggeredDeployments::SubList::deploy( Host::Human& human,
        Deployment::Method method, VaccineLimits vaccLimits )const
{
    TimeStep age = human.getAge().ts();
    if( age >= minAge && age < maxAge ){
        if( coverage >= 1.0 || util::random::bernoulli( coverage ) ){
            HumanIntervention::deploy( human, method, vaccLimits );
        }
    }
}


// ———  Derivatives of HumanInterventionComponent  ———

class RecruitmentOnlyComponent : public HumanInterventionComponent {
public:
    RecruitmentOnlyComponent( ComponentId id ) :
        HumanInterventionComponent(id, Report::MI_RECRUIT_CTS, Report::MI_RECRUIT_TIMED)
    {}
    virtual ~RecruitmentOnlyComponent() {}
    
    /// Reports to monitoring, nothing else
    virtual void deploy( Host::Human& human, Deployment::Method method, VaccineLimits ) const{
        Survey::current().addInt( reportMeasure(method), human, 1 );
    }
    
    virtual Component::Type componentType() const{
        return Component::RECRUIT_ONLY;
    }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tRecruit only";
    }
#endif
};

/// Simple treatment: no PK/PD, just remove parasites
class SimpleTreatComponent : public HumanInterventionComponent {
public:
    SimpleTreatComponent( ComponentId id, const scnXml::MDAComponent& mda ) :
        HumanInterventionComponent(id,
                Report::MI_MDA_CTS, Report::MI_MDA_TIMED)
    {
        const scnXml::Effects::OptionSequence& options = mda.getEffects().getOption();
        assert( options.size() == 1 );
        if( options[0].getPSelection() != 1.0 )
            throw util::xml_scenario_error( "sum of pSelection of a group of treatments is not 1" );
        treatId = WithinHost::WHInterface::addTreatment( options[0] );
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits ) const{
        Survey::current().addInt( reportMeasure(method), human, 1 );
        human.withinHostModel->treatment( human, treatId );
    }
    
    virtual Component::Type componentType() const{ return Component::SIMPLE_TREAT; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tTreat";
    }
#endif
    
private:
    WithinHost::TreatmentId treatId;
};

/// As SimpleTreatComponent, but with a probabilistic choice
class ProbSimpleTreatComponent : public HumanInterventionComponent {
public:
    ProbSimpleTreatComponent( ComponentId id, const scnXml::MDAComponent& mda ) :
        HumanInterventionComponent(id,
                Report::MI_MDA_CTS, Report::MI_MDA_TIMED)
    {
        const scnXml::Effects::OptionSequence& options = mda.getEffects().getOption();
        assert( options.size() >= 1 );
        treatments.reserve( options.size() );
        double cumP = 0.0;
        for( scnXml::Effects::OptionConstIterator it = options.begin(),
            end = options.end(); it != end; ++it )
        {
            cumP += it->getPSelection();
            TreatOptions treatOpts( cumP, WithinHost::WHInterface::addTreatment( *it ) );
            treatments.push_back( treatOpts );
        }
        
        // we expect the prob. to be roughly one as an error check, but allow slight deviation
        if( cumP < 0.99 || cumP > 1.01 ) throw util::xml_scenario_error( "sum of pSelection of a group of treatments is not 1" );
        for( vector<TreatOptions>::iterator it = treatments.begin(),
            end = treatments.end(); it != end; ++it )
        {
            it->cumProb /= cumP;
        }
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits ) const{
        Survey::current().addInt( reportMeasure(method), human, 1 );
        human.withinHostModel->treatment( human, selectTreatment() );
    }
    
    virtual Component::Type componentType() const{ return Component::P_S_TREAT; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tTreat";
    }
#endif
    
private:
    WithinHost::TreatmentId selectTreatment()const{
        if( treatments.size() == 1 ) return treatments[0].treatId;
        
        double x = util::random::uniform_01();      // random sample: choose
        for( vector<TreatOptions>::const_iterator it = treatments.begin(),
            end = treatments.end(); it != end; ++it )
        {
            if( it->cumProb > x ) return it->treatId;
        }
        assert( false );    // last item should have pCum=1 and x<1 in theory
        return treatments[0].treatId; // remain type safe
    }
    
    struct TreatOptions{
        double cumProb;
        WithinHost::TreatmentId treatId;
        TreatOptions( double cumProb, WithinHost::TreatmentId treatId ):
            cumProb(cumProb), treatId(treatId) {}
    };
    vector<TreatOptions> treatments;
};

HumanInterventionComponent* createSimpleTreatComponent( ComponentId id,
                                                        const scnXml::MDAComponent& mda )
{
    if( mda.getEffects().getOption().size() == 1 )
        return new SimpleTreatComponent( id, mda );
    else
        return new ProbSimpleTreatComponent( id, mda );
}

class ScreenComponent : public HumanInterventionComponent {
public:
    ScreenComponent( ComponentId id, const scnXml::Screen& elt ) :
        HumanInterventionComponent(id,
                Report::MI_SCREENING_CTS, Report::MI_SCREENING_TIMED),
        positive( elt.getPositive() ),
        negative( elt.getNegative() )
    {
        diagnostic.setXml( elt.getDiagnostic() );
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits vaccLimits ) const{
        Survey::current().addInt( reportMeasure(method), human, 1 );
        if( diagnostic.isPositive( human.withinHostModel->getTotalDensity() ) ){
            positive.deploy( human, method, vaccLimits );
        }else{
            negative.deploy( human, method, vaccLimits );
        }
    }
    
    virtual Component::Type componentType() const{ return Component::SCREEN; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tScreen";
    }
#endif
    
private:
    WithinHost::Diagnostic diagnostic;
    TriggeredDeployments positive, negative;
};

class MDA1DComponent : public HumanInterventionComponent {
public:
    MDA1DComponent( ComponentId id, const scnXml::DecisionTree& description ) :
        HumanInterventionComponent(id,
                Report::MI_MDA_CTS, Report::MI_MDA_TIMED),
        m_screenMeasureCts(Report::MI_SCREENING_CTS),
        m_screenMeasureTimed(Report::MI_SCREENING_TIMED)
    {
        if( !util::ModelOptions::option( util::CLINICAL_EVENT_SCHEDULER ) )
            throw util::xml_scenario_error( "MDA1D intervention: requires CLINICAL_EVENT_SCHEDULER option" );
        Clinical::ESCaseManagement::initMDA( description );
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits ) const{
        human.getClinicalModel().massDrugAdministration( human, screeningMeasure(method), reportMeasure(method) );
    }
    
    virtual Component::Type componentType() const{ return Component::MDA_TS1D; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tMDA1D";
    }
#endif
    
private:
    /** Trivial helper function to get deployment measure. */
    inline ReportMeasureI screeningMeasure( Deployment::Method method )const{
        return (method == Deployment::TIMED) ? m_screenMeasureTimed :
            (method == Deployment::CTS) ? m_screenMeasureCts :
            Report::MI_TREAT_DEPLOYMENTS;
    }
    ReportMeasureI m_screenMeasureCts, m_screenMeasureTimed;
};

class ClearImmunityComponent : public HumanInterventionComponent {
public:
    ClearImmunityComponent( ComponentId id ) :
        HumanInterventionComponent(id, Report::MI_NUM, Report::MI_NUM /*never reported*/) {}
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits )const{
        human.clearImmunity();
    }
    
    virtual Component::Type componentType() const{ return Component::CLEAR_IMMUNITY; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tclear immunity";
    }
#endif
};

} }
