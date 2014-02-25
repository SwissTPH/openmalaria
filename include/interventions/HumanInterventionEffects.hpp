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
#include "WithinHost/DescriptiveIPTWithinHost.h"
#include <schema/healthSystem.h>
#include <schema/interventions.h>

namespace OM { namespace interventions {
    using Host::Human;

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
    MDAEffect( EffectId id, const scnXml::MDA& mda ) : HumanInterventionEffect(id) {
        if( !mda.getDiagnostic().present() ){
            // Note: allow no description for now to avoid XML changes.
            //throw util::xml_scenario_error( "error: interventions.MDA.diagnostic element required for MDA with 5-day timestep" );
            diagnostic.setDeterministic( 0.0 );
        }else{
            diagnostic.setXml( mda.getDiagnostic().get() );
        }
        const scnXml::DrugWithCompliance& drug = mda.getDrugEffect();
        double pCompliance = drug.getCompliance().getPCompliance();
        double nonComplierMult = drug.getCompliance().getNonCompliersMultiplier();
        double mult = pCompliance + (1.0 - pCompliance) * nonComplierMult;
        const scnXml::CompliersEffective::TimestepSequence& seq =
                drug.getCompliersEffective().getTimestep();
        size_t len = seq.size();
        pClearanceByTime.resize(len);
        for( size_t i = 0; i < len; ++i ){
            double pCompliers = seq[i].getPClearance();
            pClearanceByTime[i] = mult * pCompliers;
        }
        if( len < 1 ){
            throw util::xml_scenario_error(
                "interventions.human.effect.MDA.drugEffect: require at "
                "least one timestep element" );
        }
        if( len > 1 ){
            if( true /*model not implemented*/ )
                throw util::xml_scenario_error(
                    "MDA with prophylactic effect (unimplemented)" );
        }   
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits ) const{
        //TODO(monitoring): separate reports for mass and continuous deployments
        
        Monitoring::Surveys.getSurvey(human.isInAnyCohort())
                .reportMassScreening(human.getMonitoringAgeGroup(), 1);
        if( !diagnostic.isPositive( human.withinHostModel->getTotalDensity() ) ){
            return;
        }
        Monitoring::Surveys.getSurvey(human.isInAnyCohort())
                .reportMDA(human.getMonitoringAgeGroup(), 1);
	
        double pClearance = pClearanceByTime[0];
        if( pClearance >= 1.0 || util::random::bernoulli( pClearance ) ){
            human.withinHostModel->clearInfections(false/*value doesn't matter*/);
        }
        if( pClearanceByTime.size() > 1 ){
            //TODO:...
//             human.withinHostModel->addProphylacticEffects( pClearanceByTime );
        }
    }
    
    virtual Effect::Type effectType() const{ return Effect::MDA; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tMDA";
    }
#endif
    
private:
    WithinHost::Diagnostic diagnostic;
    vector<double> pClearanceByTime;
};

class MDA1DEffect : public HumanInterventionEffect {
public:
    MDA1DEffect( EffectId id, const scnXml::MDA1D& description ) : HumanInterventionEffect(id) {
	if( !util::ModelOptions::option( util::CLINICAL_EVENT_SCHEDULER ) )
	  throw util::xml_scenario_error( "MDA1D intervention: requires CLINICAL_EVENT_SCHEDULER option" );
        Clinical::ESCaseManagement::initMDA( description );
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits ) const{
        human.getClinicalModel().massDrugAdministration ( human );
    }
    
    virtual Effect::Type effectType() const{ return Effect::MDA_TS1D; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\tMDA1D";
    }
#endif
    
private:
};

class VaccineEffect : public HumanInterventionEffect {
public:
    VaccineEffect( EffectId id, const scnXml::VaccineDescription& seq, Vaccine::Types type ) :
            HumanInterventionEffect(id), type(type)
    {
        Vaccine::init( seq, type );
    }
    
    void deploy( Human& human, Deployment::Method method, VaccineLimits vaccLimits )const{
        human.getVaccine().possiblyVaccinate( human, method, type, vaccLimits );
    }
    
    virtual Effect::Type effectType() const{
        if( type == Vaccine::PEV ) return Effect::PEV;
        else if( type == Vaccine::BSV ) return Effect::BSV;
        else if( type == Vaccine::TBV ) return Effect::TBV;
        else throw SWITCH_DEFAULT_EXCEPTION;
    }
    
    inline Vaccine::Types getVaccineType()const{ return type; }
    
#ifdef WITHOUT_BOINC
    virtual void print_details( std::ostream& out )const{
        out << id().id << "\t" <<
            (type == Vaccine::PEV ? "PEV" : (type == Vaccine::BSV ? "BSV" : "TBV"))
           ;
    }
#endif
    
private:
    Vaccine::Types type;
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
