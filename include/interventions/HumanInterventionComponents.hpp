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
#include "interventions/Cohort.h"
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

// ———  HumanInterventionComponent  ———

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

bool componentCmp(const HumanInterventionComponent *a, const HumanInterventionComponent *b){
    return a->componentType() < b->componentType();
}
void HumanIntervention::sortComponents(){
    std::stable_sort( components.begin(), components.end(), componentCmp );
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

// ———  Derivatives  ———

class MDAComponentBase : public HumanInterventionComponent {
protected:
    MDAComponentBase( ComponentId id ) :
        HumanInterventionComponent(id,
                Report::MI_MDA_CTS, Report::MI_MDA_TIMED),
        m_screenMeasureCts(Report::MI_SCREENING_CTS),
        m_screenMeasureTimed(Report::MI_SCREENING_TIMED)
    {}
    
    /** Trivial helper function to get deployment measure. */
    inline ReportMeasureI screeningMeasure( Deployment::Method method )const{
        return (method == Deployment::TIMED) ? m_screenMeasureTimed : m_screenMeasureCts;
    }
    
private:
    ReportMeasureI m_screenMeasureCts, m_screenMeasureTimed;
};

class MDAComponent : public MDAComponentBase {
public:
    MDAComponent( ComponentId id, const scnXml::MDAComponent& mda ) :
        MDAComponentBase(id)
    {
        if( !mda.getDiagnostic().present() ){
            // Note: allow no description for now to avoid XML changes.
            //throw util::xml_scenario_error( "error: interventions.MDA.diagnostic element required for MDA with 5-day timestep" );
            diagnostic.setDeterministic( 0.0 );
        }else{
            diagnostic.setXml( mda.getDiagnostic().get() );
        }
        
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
        //TODO(monitoring): separate reports for mass and continuous deployments
        
        Survey::current().addInt( screeningMeasure(method), human, 1 );
        if( !diagnostic.isPositive( human.withinHostModel->getTotalDensity() ) ){
            return;
        }
        Survey::current().addInt( reportMeasure(method), human, 1 );
        
        human.withinHostModel->treatment( selectTreatment() );
    }
    
    virtual Component::Type componentType() const{ return Component::MDA; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tMDA";
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
    
    WithinHost::Diagnostic diagnostic;
    struct TreatOptions{
        double cumProb;
        WithinHost::TreatmentId treatId;
        TreatOptions( double cumProb, WithinHost::TreatmentId treatId ):
            cumProb(cumProb), treatId(treatId) {}
    };
    vector<TreatOptions> treatments;
};

class MDA1DComponent : public MDAComponentBase {
public:
    MDA1DComponent( ComponentId id, const scnXml::HSESCaseManagement& description ) :
        MDAComponentBase(id)
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
