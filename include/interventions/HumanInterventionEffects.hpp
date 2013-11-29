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

// Note: this file is included by exactly one source (Interventions.cpp)
// and contains definitions as well as declarations.

// The includes here are more for documentation than required.
#include "interventions/Interventions.h"
#include "interventions/Cohort.h"
#include "Monitoring/Surveys.h"
#include "util/ModelOptions.h"
#include "Clinical/ESCaseManagement.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Clinical/Diagnostic.h"
#include "Clinical/CaseManagementCommon.h"
#include "Host/Human.h"
#include "Transmission/TransmissionModel.h"
#include "WithinHost/DescriptiveIPTWithinHost.h"
#include <schema/healthSystem.h>
#include <schema/interventions.h>

namespace OM { namespace interventions {
    using Host::Human;

class MDAEffect : public HumanInterventionEffect {
public:
    MDAEffect( size_t index, const scnXml::MDA& mda ) : HumanInterventionEffect(index) {
        if( !mda.getDiagnostic().present() ){
            // Note: allow no description for now to avoid XML changes.
            //throw util::xml_scenario_error( "error: interventions.MDA.diagnostic element required for MDA with 5-day timestep" );
            scnXml::HSDiagnostic diagnostic;
            scnXml::Deterministic det(0.0);
            diagnostic.setDeterministic(det);
            this->diagnostic.init(diagnostic);
        }else{
            diagnostic.init( mda.getDiagnostic().get() );
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
            if( !util::ModelOptions::option( util::PROPHYLACTIC_DRUG_ACTION_MODEL ) )
                throw util::xml_scenario_error(
                    "MDA with prophylactic effect (i.e. with more than one"
                    " timestep element in drugEffect element) requires the"
                    " PROPHYLACTIC_DRUG_ACTION_MODEL" );
        }   
    }
    
    void deploy( Human& human, Deployment::Method method ) const{
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
            human.withinHostModel->clearAllInfections();
        }
        if( pClearanceByTime.size() > 1 )
            human.withinHostModel->addProphylacticEffects( pClearanceByTime );
    }
    
    virtual Effect::Type effectType() const{ return Effect::MDA; }
    
private:
    Clinical::Diagnostic diagnostic;
    vector<double> pClearanceByTime;
};

class MDA1DEffect : public HumanInterventionEffect {
public:
    MDA1DEffect( size_t index, const scnXml::MDA1D& description ) : HumanInterventionEffect(index) {
	if( !util::ModelOptions::option( util::CLINICAL_EVENT_SCHEDULER ) )
	  throw util::xml_scenario_error( "MDA1D intervention: requires CLINICAL_EVENT_SCHEDULER option" );
        Clinical::ESCaseManagement::initMDA( description );
    }
    
    void deploy( Human& human, Deployment::Method method ) const{
        human.getClinicalModel().massDrugAdministration ( human );
    }
    
    virtual Effect::Type effectType() const{ return Effect::MDA_TS1D; }
    
private:
};

class VaccineEffect : public HumanInterventionEffect {
public:
    VaccineEffect( size_t index, const scnXml::VaccineDescription& seq, Host::Vaccine::Types type ) :
            HumanInterventionEffect(index), type(type)
    {
        Host::Vaccine::init( seq, type );
    }
    
    void deploy( Human& human, Deployment::Method method )const{
        human.getVaccine().vaccinate( human, method, type );
    }
    
    virtual Effect::Type effectType() const{
        if( type == Host::Vaccine::PEV ) return Effect::PEV;
        else if( type == Host::Vaccine::BSV ) return Effect::BSV;
        else if( type == Host::Vaccine::TBV ) return Effect::TBV;
        else throw SWITCH_DEFAULT_EXCEPTION;
    }
    
    inline Host::Vaccine::Types getVaccineType()const{ return type; }
    
private:
    Host::Vaccine::Types type;
};

class ClearImmunityEffect : public HumanInterventionEffect {
public:
    ClearImmunityEffect( size_t index ) : HumanInterventionEffect(index) {}
    
    void deploy( Human& human, Deployment::Method method )const{
        human.clearImmunity();
    }
    
    virtual Effect::Type effectType() const{ return Effect::CLEAR_IMMUNITY; }
};

} }
