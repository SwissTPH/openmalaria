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
#include "Monitoring/Surveys.h"
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
    using Monitoring::Survey;

// ———  HumanInterventionEffect  ———

void HumanIntervention::deploy( Human& human, Deployment::Method method,
    VaccineLimits vaccLimits ) const
{
    for( vector<const HumanInterventionEffect*>::const_iterator it = effects.begin();
            it != effects.end(); ++it )
    {
        const interventions::HumanInterventionEffect& effect = **it;
        effect.deploy( human, method, vaccLimits );
        human.updateLastDeployed( effect.id() );
    }
}

bool effectCmp(const HumanInterventionEffect *a, const HumanInterventionEffect *b){
    return a->effectType() < b->effectType();
}
void HumanIntervention::sortEffects(){
    std::stable_sort( effects.begin(), effects.end(), effectCmp );
}

#ifdef WITHOUT_BOINC
void HumanIntervention::print_details( std::ostream& out )const{
    out << "human:";
    for( vector<const HumanInterventionEffect*>::const_iterator it =
        effects.begin(); it != effects.end(); ++it ){
        out << '\t' << (*it)->id().id;
    }
}
#endif

// ———  Derivatives  ———

class MDAEffect : public HumanInterventionEffect {
public:
    MDAEffect( EffectId id, const scnXml::MDA& mda ) :
        HumanInterventionEffect(id)
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
        
        Survey& survey = Monitoring::Surveys.getSurvey(human.isInAnyCohort());
        survey.addInt( (method == Deployment::TIMED) ? Survey::MI_SCREENING_TIMED :
                       Survey::MI_SCREENING_CTS, human.getMonitoringAgeGroup(), 1 );
        if( !diagnostic.isPositive( human.withinHostModel->getTotalDensity() ) ){
            return;
        }
        survey.addInt( (method == Deployment::TIMED) ? Survey::MI_MDA_TIMED :
                       Survey::MI_MDA_CTS, human.getMonitoringAgeGroup(), 1 );
        
        human.withinHostModel->treatment( selectTreatment() );
    }
    
    virtual Effect::Type effectType() const{ return Effect::MDA; }
    
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

class MDA1DEffect : public HumanInterventionEffect {
public:
    MDA1DEffect( EffectId id, const scnXml::HSESCaseManagement& description ) : HumanInterventionEffect(id) {
	if( !util::ModelOptions::option( util::CLINICAL_EVENT_SCHEDULER ) )
	  throw util::xml_scenario_error( "MDA1D intervention: requires CLINICAL_EVENT_SCHEDULER option" );
        Clinical::ESCaseManagement::initMDA( description );
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits ) const{
        human.getClinicalModel().massDrugAdministration( method, human );
    }
    
    virtual Effect::Type effectType() const{ return Effect::MDA_TS1D; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tMDA1D";
    }
#endif
    
private:
};

class ClearImmunityEffect : public HumanInterventionEffect {
public:
    ClearImmunityEffect( EffectId id ) : HumanInterventionEffect(id) {}
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits )const{
        human.clearImmunity();
    }
    
    virtual Effect::Type effectType() const{ return Effect::CLEAR_IMMUNITY; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tclear immunity";
    }
#endif
};

} }
