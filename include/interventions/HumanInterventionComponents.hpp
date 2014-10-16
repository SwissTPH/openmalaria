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
#include "PkPd/LSTMTreatments.h"
#include "util/random.h"
#include <schema/healthSystem.h>
#include <schema/interventions.h>

namespace OM { namespace interventions {
using Host::Human;
using WithinHost::Diagnostic;
using WithinHost::diagnostics;
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
        human.reportDeployment( component.id(), component.duration() );
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
        minAge( sim::fromYearsN( elt.getMinAge() ) ),
        maxAge( sim::future() ),
        coverage( elt.getP() )
{
    if( elt.getMaxAge().present() )
        maxAge = sim::fromYearsN( elt.getMaxAge().get() );
    
    if( minAge < sim::zero() || maxAge < minAge ){
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
    // This may be used within a time step; in that case we use human age at
    // the beginning of the step since this gives the expected result with age
    // limits of 0 to 1 time step.
    SimTime age = human.age(sim::nowOrTs1()/*TODO: should be now() but requires delayed triggered deployments*/);
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

class ScreenComponent : public HumanInterventionComponent {
public:
    ScreenComponent( ComponentId id, const scnXml::Screen& elt ) :
        HumanInterventionComponent(id,
                Report::MI_SCREENING_CTS, Report::MI_SCREENING_TIMED),
        diagnostic(diagnostics::get(elt.getDiagnostic())),
        positive( elt.getPositive() ),
        negative( elt.getNegative() )
    {}
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits vaccLimits ) const{
        Survey::current().addInt( reportMeasure(method), human, 1 );
        if( human.withinHostModel->diagnosticResult(diagnostic) ){
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
    const Diagnostic& diagnostic;
    HumanIntervention positive, negative;
};

/// Simple treatment: no PK/PD, just remove parasites
class TreatSimpleComponent : public HumanInterventionComponent {
public:
    TreatSimpleComponent( ComponentId id, const scnXml::DTTreatSimple& elt ) :
        HumanInterventionComponent(id,
                Report::MI_MDA_CTS, Report::MI_MDA_TIMED)
    {
        //NOTE: this code is currently identical to that in CMDTTreatSimple
        try{
            SimTime durL = UnitParse::readShortDuration( elt.getDurationLiver(), UnitParse::NONE ),
                durB = UnitParse::readShortDuration( elt.getDurationBlood(), UnitParse::NONE );
            SimTime neg1 = -sim::oneTS();
            if( durL < neg1 || durB < neg1 ){
                throw util::xml_scenario_error( "treatSimple: cannot have durationBlood or durationLiver less than -1" );
            }
            if( util::ModelOptions::option( util::VIVAX_SIMPLE_MODEL ) ){
                if( durL != sim::zero() || durB != neg1 )
                    throw util::unimplemented_exception( "vivax model only supports timestepsLiver=0, timestepsBlood=-1" );
                // Actually, the model ignores these parameters; we just don't want somebody thinking it doesn't.
            }
            timeLiver = durL;
            timeBlood = durB;
        }catch( const util::format_error& e ){
            throw util::xml_scenario_error( string("treatSimple: ").append(e.message()) );
        }
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits ) const{
        Survey::current().addInt( reportMeasure(method), human, 1 );
        human.withinHostModel->treatSimple( timeLiver, timeBlood );
    }
    
    virtual Component::Type componentType() const{ return Component::TREAT_SIMPLE; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\ttreatSimple";
    }
#endif
    
private:
    SimTime timeLiver, timeBlood;
};

//NOTE: this is similar to CMDTTreatPKPD, but supporting only a single "treatment"
class TreatPKPDComponent : public HumanInterventionComponent {
public:
    TreatPKPDComponent( ComponentId id, const scnXml::DTTreatPKPD& elt ) :
        HumanInterventionComponent(id,
                Report::MI_MDA_CTS, Report::MI_MDA_TIMED),
            schedule(PkPd::LSTMTreatments::findSchedule(elt.getSchedule())),
            dosage(PkPd::LSTMTreatments::findDosages(elt.getDosage())),
            delay_h(elt.getDelay_h())
    {
        if( !util::ModelOptions::option( util::INCLUDES_PK_PD ) ){
            throw util::xml_scenario_error( "treatPKPD: requires INCLUDES_PK_PD model option" );
        }
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits ) const{
        Survey::current().addInt( reportMeasure(method), human, 1 );
        human.withinHostModel->treatPkPd( schedule, dosage, delay_h );
    }
    
    virtual Component::Type componentType() const{ return Component::TREAT_PKPD; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\ttreatPKPD";
    }
#endif
    
private:
    size_t schedule;        // index of the schedule
    size_t dosage;          // name of the dosage table
    double delay_h;         // delay in hours
};

class DecisionTreeComponent : public HumanInterventionComponent {
public:
    DecisionTreeComponent( ComponentId id, const scnXml::DecisionTree& elt ) :
        HumanInterventionComponent(id, Report::MI_MDA_CTS, Report::MI_MDA_TIMED),
        tree(Clinical::CMDecisionTree::create(elt, false))
    {}
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits )const{
        Clinical::CMDTOut out = tree.exec( Clinical::CMHostData(human,
                human.age(sim::nowOrTs1()).inYears(),
                Clinical::Episode::NONE /*parameter not needed*/) );
        if( out.treated ){
            Survey::current().addInt( reportMeasure(method), human, 1 );
        }
    }
    
    virtual Component::Type componentType() const{ return Component::CM_DT; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tDecision tree";
    }
#endif
    
private:
    const Clinical::CMDecisionTree& tree;
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
