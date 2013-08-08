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

#include "interventions/Interventions.h"
#include "Host/Vaccine.h"
#include "Population.h"
#include "util/random.h"
#include "Clinical/ESCaseManagement.h"
#include "Clinical/ImmediateOutcomes.h"
#include "Clinical/Diagnostic.h"
#include "WithinHost/DescriptiveIPTWithinHost.h"
#include "Clinical/CaseManagementCommon.h"
#include "Monitoring/Surveys.h"
#include "interventions/GVI.h"

namespace OM { namespace interventions {
    using Host::Human;

// ———  ContinuousDeployment and derivatives  ———

ContinuousDeployment::ContinuousDeployment(
        const ::scnXml::ContinuousDeployment& elt ) :
    begin( elt.getBegin() ),
    end( elt.getEnd() ),
    deployAge( TimeStep::fromYears( elt.getTargetAgeYrs() ) ),
    cohortOnly( elt.getCohort() ),
    coverage( elt.getCoverage() )
{
    if( begin < TimeStep(0) || end < begin ){
        throw util::xml_scenario_error("continuous intervention must have 0 <= begin <= end");
    }
    if( deployAge <= TimeStep(0) ){
        ostringstream msg;
        msg << "continuous intervention with target age "<<elt.getTargetAgeYrs();
        msg << " years corresponds to timestep "<<deployAge;
        msg << "; must be at least timestep 1.";
        throw util::xml_scenario_error( msg.str() );
    }
    if( deployAge > TimeStep::maxAgeIntervals ){
        ostringstream msg;
        msg << "continuous intervention must have target age no greater than ";
        msg << TimeStep::maxAgeIntervals * TimeStep::yearsPerInterval;
        throw util::xml_scenario_error( msg.str() );
    }
    if( !(coverage >= 0.0 && coverage <= 1.0) ){
        throw util::xml_scenario_error("continuous intervention coverage must be in range [0,1]");
    }
}

bool ContinuousDeployment::filterAndDeploy( Host::Human& human, const Population& population )const{
    TimeStep age = TimeStep::simulation - human.getDateOfBirth();
    if( deployAge > age ){
        // stop processing continuous deployments for this
        // human for now because remaining ones happen in the future
        return false;
    }else if( deployAge == age ){
        if( begin <= TimeStep::interventionPeriod &&
            TimeStep::interventionPeriod < end &&
            ( !cohortOnly || human.isInCohort() ) &&
            util::random::uniform_01() < coverage )     // RNG call should be last test
        {
            deploy( human, population );
        }
    }//else: for some reason, a deployment age was missed; ignore it
    return true;
}

/** Age-based deployment for new list-of-effect interventions. */
class ContinuousHumanIntervention : public ContinuousDeployment {
public:
    ContinuousHumanIntervention( const scnXml::ContinuousDeployment& elt,
                                 const HumanIntervention* intervention ) :
        ContinuousDeployment( elt ),
        intervention( intervention )
    {
    }
    
protected:
    virtual void deploy( Host::Human& human, const Population& population )const{
        intervention->deploy( human, Deployment::CTS );
    }
    
    const HumanIntervention *intervention;
};


// ———  TimedDeployment and derivatives  ———

TimedDeployment::TimedDeployment(TimeStep deploymentTime) :
    time( deploymentTime )
{
    if( deploymentTime < TimeStep(0) ){
        throw util::xml_scenario_error("timed intervention deployment: may not be negative");
    }else if( deploymentTime >= Monitoring::Surveys.getFinalTimestep() ){
        cerr << "Warning: timed intervention deployment at time "<<deploymentTime.asInt();
        cerr << " happens after last survey" << endl;
    }
}

class DummyTimedDeployment : public TimedDeployment {
public:
    DummyTimedDeployment() :
        TimedDeployment( TimeStep(0) )
    {
        // TimedDeployment's ctor checks that the deployment time-step is
        // within the intervention period. We want this time to be after the
        // last time-step, so set the time here after TimedDeployment's ctor
        // check has been done (hacky).
        time = TimeStep::future;
    }
    virtual void deploy (OM::Population&) {}
};

class TimedChangeHSDeployment : public TimedDeployment {
public:
    TimedChangeHSDeployment( const scnXml::ChangeHS::TimedDeploymentType& hs ) :
        TimedDeployment( TimeStep( hs.getTime() ) ),
        newHS( hs._clone() )
    {}
    virtual void deploy (OM::Population& population) {
        Clinical::CaseManagementCommon::changeHealthSystem( *newHS );
        delete newHS;
        newHS = 0;
    }
    
private:
    scnXml::HealthSystem *newHS;
};

class TimedChangeEIRDeployment : public TimedDeployment {
public:
    TimedChangeEIRDeployment( const scnXml::ChangeEIR::TimedDeploymentType& nv ) :
        TimedDeployment( TimeStep( nv.getTime() ) ),
        newEIR( nv._clone() )
    {}
    virtual void deploy (OM::Population& population) {
        population.transmissionModel().changeEIRIntervention( *newEIR );
        delete newEIR;
        newEIR = 0;
    }
    
private:
    scnXml::NonVector *newEIR;
};

class TimedUninfectVectorsDeployment : public TimedDeployment {
public:
    TimedUninfectVectorsDeployment( TimeStep deployTime ) :
        TimedDeployment( deployTime )
    {}
    virtual void deploy (OM::Population& population) {
        population.transmissionModel().uninfectVectors();
    }
};

class TimedR_0Deployment : public TimedDeployment {
public:
    TimedR_0Deployment( TimeStep deployTime ) :
        TimedDeployment( deployTime )
    {}
    virtual void deploy (OM::Population& population) {
        int i = (int)std::floor (util::random::uniform_01() * population.size());        // pick a human
        Population::Iter it = population.begin();
        while (i > 0){  // find human (can't use population[i])
            ++it;
            --i;
        }
        assert( i == 0 );
        assert( it != population.end() );
        it->R_0Vaccines();
        it->addInfection();
    }
};

/// Timed deployment of human-specific interventions
class TimedHumanDeployment : public TimedDeployment {
public:
    /** 
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy. */
    TimedHumanDeployment( const scnXml::Mass& mass,
                           const HumanIntervention* intervention ) :
        TimedDeployment( TimeStep( mass.getTime() ) ),
        minAge( TimeStep::fromYears( mass.getMinAge() ) ),
        maxAge( TimeStep::fromYears( mass.getMaxAge() ) ),
        cohortOnly( mass.getCohort() ),
        coverage( mass.getCoverage() ),
        intervention( intervention )
    {
        if( !(coverage >= 0.0 && coverage <= 1.0) ){
            throw util::xml_scenario_error("timed intervention coverage must be in range [0,1]");
        }
        if( minAge < TimeStep(0) || maxAge < minAge ){
            throw util::xml_scenario_error("timed intervention must have 0 <= minAge <= maxAge");
        }
    }
    
    virtual void deploy (OM::Population& population) {
        for (Population::Iter iter = population.begin(); iter != population.end(); ++iter) {
            TimeStep age = TimeStep::simulation - iter->getDateOfBirth();
            if( age >= minAge && age < maxAge ){
                if( !cohortOnly || iter->isInCohort() ){
                    if( util::random::uniform_01() < coverage ){
                        intervention->deploy( *iter, Deployment::TIMED );
                    }
                }
            }
        }
    }
    
protected:
    // restrictions on deployment
    TimeStep minAge;
    TimeStep maxAge;
    bool cohortOnly;
    double coverage;    // proportion coverage within group meeting above restrictions
    const HumanIntervention *intervention;
};

/// Timed deployment of human-specific interventions in cumulative mode
class TimedCumulativeHumanDeployment : public TimedHumanDeployment {
public:
    /** 
     * @param mass XML element specifying the age range and compliance
     * (proportion of eligible individuals who receive the intervention).
     * @param intervention The HumanIntervention to deploy.
     * @param effect_index Index of effect to test coverage for
     * @param maxAge Maximum time-span to consider a deployed effect still to be effective */
    TimedCumulativeHumanDeployment( const scnXml::Mass& mass,
                           const HumanIntervention* intervention,
                           size_t effect_index, TimeStep maxAge ) :
        TimedHumanDeployment( mass, intervention ),
        cumCovInd( effect_index ), maxInterventionAge( maxAge )
    {
    }
    
    virtual void deploy (OM::Population& population) {
        // Cumulative case: bring target group's coverage up to target coverage
        vector<Host::Human*> unprotected;
        size_t total = 0;       // number of humans within age bound and optionally cohort
        for (Population::Iter iter = population.begin(); iter != population.end(); ++iter) {
            TimeStep age = TimeStep::simulation - iter->getDateOfBirth();
            if( age >= minAge && age < maxAge ){
                if( !cohortOnly || iter->isInCohort() ){
                    total+=1;
                    if( iter->needsRedeployment(cumCovInd, maxInterventionAge) )
                        unprotected.push_back( &*iter );
                }
            }
        }
        
        double propProtected = static_cast<double>( total - unprotected.size() ) / static_cast<double>( total );
        if( propProtected < coverage ){
            // Proportion propProtected are already covered, so need to
            // additionally cover the proportion (coverage - propProtected),
            // selected from the list unprotected.
            double additionalCoverage = (coverage - propProtected) / (1.0 - propProtected);
            for (vector<Host::Human*>::iterator iter = unprotected.begin();
                 iter != unprotected.end(); ++iter)
            {
                if( util::random::uniform_01() < additionalCoverage ){
                    intervention->deploy( **iter, Deployment::TIMED );
                }
            }
        }
    }
    
protected:
    size_t cumCovInd;
    // max age at which an intervention is considered not to need replacement
    TimeStep maxInterventionAge;
};

class TimedVectorDeployment : public TimedDeployment {
public:
    TimedVectorDeployment( TimeStep deployTime, size_t instance ) :
        TimedDeployment( deployTime ),
        inst(instance)
    {}
    virtual void deploy (OM::Population& population) {
      population.transmissionModel().deployVectorPopInterv(inst);
    }
private:
    size_t inst;
};

// ———  HumanInterventionEffect  ———

void HumanIntervention::deploy( Human& human, Deployment::Method method ) const{
    for( vector<const HumanInterventionEffect*>::const_iterator it = effects.begin();
            it != effects.end(); ++it )
    {
        human.deploy( **it, method );
    }
}

bool effectCmp(const HumanInterventionEffect *a, const HumanInterventionEffect *b){
    return a->effectType() < b->effectType();
}
void HumanIntervention::sortEffects(){
    std::stable_sort( effects.begin(), effects.end(), effectCmp );
}

class MDAEffect : public HumanInterventionEffect {
public:
    MDAEffect( size_t index, const scnXml::MDA& mda ) : HumanInterventionEffect(index) {
        if( TimeStep::interval == 5 ){
            if( !mda.getDrugEffect().present() )
                throw util::xml_scenario_error( "interventions.human.effect.MDA: drugEffect element required for 5-day timestep" );
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
            const scnXml::DrugWithCompliance& drug = mda.getDrugEffect().get();
            double pCompliance = drug.getCompliance().getPCompliance();
            double nonComplierMult = drug.getCompliance().getNonCompliersMultiplier();
            double mult = pCompliance + (1.0 - pCompliance) * nonComplierMult;
            const scnXml::CompliersEffective::TimestepSequence& seq = drug.getCompliersEffective().getTimestep();
            size_t len = seq.size();
            pClearanceByTime.resize(len);
            for( size_t i = 0; i < len; ++i ){
                double pCompliers = seq[i].getPClearance();
                pClearanceByTime[i] = mult * pCompliers;
            }
            if( len < 1 ){
                throw util::xml_scenario_error( "interventions.human.effect.MDA.drugEffect: require at least one timestep element" );
            }
            if( len > 1 ){
                throw util::unimplemented_exception( "MDA with prophylactic effect (i.e. more than one timestep element in drugEffect element" );
            }
        }else{
            // We could either do what we did before for the 1-day timestep
            // (below, and calling ClinicalEventScheduler::massDrugAdministration
            // on deployment), or we could use the 5-day timestep approach
            // (or even support both).
            throw util::unimplemented_exception( "MDA/MSAT on 1-day timestep (disabled pending review)" );
            if( !mda.getDescription().present() ){
                throw util::xml_scenario_error( "error: interventions.MDA.description element required for MDA with 1-day timestep" );
            }
            Clinical::ESCaseManagement::initMDA( mda.getDescription().get() );
        }
    }
    
    void deploy( Human& human, Deployment::Method method ) const{
        //TODO: shouldn't really use the same reports for mass and continuous deployment, right?
        
        Monitoring::Surveys.getSurvey(human.isInCohort()).reportMassScreening(human.getMonitoringAgeGroup(), 1);
        if( !diagnostic.isPositive( human.withinHostModel->getTotalDensity() ) ){
            return;
        }
        double pClearance = pClearanceByTime[0];
        if( pClearance >= 1.0 || util::random::bernoulli( pClearance ) ){
            human.withinHostModel->clearInfections(human.getClinicalModel().latestIsSevere());
        }
        Monitoring::Surveys.getSurvey(human.isInCohort()).reportMDA(human.getMonitoringAgeGroup(), 1);
    }
    
    virtual Effect::Type effectType() const{ return Effect::MDA; }
    
private:
    Clinical::Diagnostic diagnostic;
    vector<double> pClearanceByTime;
};

class VaccineEffect : public HumanInterventionEffect {
public:
    VaccineEffect( size_t index, const scnXml::VaccineDescription& seq, Host::Vaccine::Types type ) :
            HumanInterventionEffect(index), type(type)
    {
        Host::Vaccine::types[type].initVaccine( seq, type );
    }
    
    void deploy( Human& human, Deployment::Method method )const{
        human.deployVaccine( method, type );
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

class IPTEffect : public HumanInterventionEffect {
public:
    IPTEffect( size_t index, const scnXml::IPTDescription& elt ) : HumanInterventionEffect(index){
        // Check compatibilities with MDA model, in particular use of clearInfections(bool).
        // Is it needed now MDA allows prophylactic effects and continuous deployment, anyway?
        throw util::unimplemented_exception( "IPT model (disabled pending review)" );
        WithinHost::DescriptiveIPTWithinHost::init( elt );
    }
    
    void deploy( Human& human, Deployment::Method method )const{
        human.deployIPT( method );
    }
    
    virtual Effect::Type effectType() const{ return Effect::IPT; }
};

class ITNEffect : public HumanInterventionEffect {
public:
    ITNEffect( size_t index, const scnXml::ITNDescription& elt,
               Transmission::TransmissionModel& transmissionModel ) : HumanInterventionEffect(index),
               transmission( transmissionModel )
    {
        transmissionModel.setITNDescription( elt );
    }
    
    void deploy( Human& human, Deployment::Method method )const{
        human.deployITN( method, transmission );
    }
    
    virtual Effect::Type effectType() const{ return Effect::ITN; }
    
private:
    Transmission::TransmissionModel& transmission;      //TODO: storing this is not a nice solution; do we need to pass?
};

class IRSEffect : public HumanInterventionEffect {
public:
    IRSEffect( size_t index, const scnXml::IRSDescription& elt,
               Transmission::TransmissionModel& transmissionModel ) : HumanInterventionEffect(index),
               transmission( transmissionModel )
    {
        transmissionModel.setIRSDescription( elt );
    }
    
    void deploy( Human& human, Deployment::Method method )const{
        human.deployIRS( method, transmission );
    }
    
    virtual Effect::Type effectType() const{ return Effect::IRS; }
    
private:
    Transmission::TransmissionModel& transmission;      //TODO: storing this is not a nice solution; do we need to pass?
};

class CohortSelectionEffect : public HumanInterventionEffect {
public:
    CohortSelectionEffect( size_t index ) : HumanInterventionEffect(index) {}
    
    void deploy( Human& human, Deployment::Method method )const{
        human.addToCohort();
    }
    
    virtual Effect::Type effectType() const{ return Effect::COHORT; }
};

class ClearImmunityEffect : public HumanInterventionEffect {
public:
    ClearImmunityEffect( size_t index ) : HumanInterventionEffect(index) {}
    
    void deploy( Human& human, Deployment::Method method )const{
        human.clearImmunity();
    }
    
    virtual Effect::Type effectType() const{ return Effect::CLEAR_IMMUNITY; }
};


// ———  InterventionManager  ———

InterventionManager::InterventionManager (const scnXml::Interventions& intervElt, OM::Population& population) :
    nextTimed(0), _cohortEnabled(false)
{
    if( intervElt.getChangeHS().present() ){
        const scnXml::ChangeHS& chs = intervElt.getChangeHS().get();
        if( chs.getTimedDeployment().size() > 0 ){
            // timed deployments:
            typedef scnXml::ChangeHS::TimedDeploymentSequence::const_iterator It;
            for( It it = chs.getTimedDeployment().begin(); it != chs.getTimedDeployment().end(); ++it ){
                timed.push_back( new TimedChangeHSDeployment( *it ) );
            }
        }
    }
    if( intervElt.getChangeEIR().present() ){
        const scnXml::ChangeEIR& eir = intervElt.getChangeEIR().get();
        if( eir.getTimedDeployment().size() > 0 ){
            // timed deployments:
            typedef scnXml::ChangeEIR::TimedDeploymentSequence::const_iterator It;
            for( It it = eir.getTimedDeployment().begin(); it != eir.getTimedDeployment().end(); ++it ){
                timed.push_back( new TimedChangeEIRDeployment( *it ) );
            }
        }
    }
    Transmission::TransmissionModel& transmission = population.transmissionModel();
    const map<string,size_t>* species_index_map = 0;
    if( intervElt.getHuman().present() ){
        const scnXml::HumanInterventions& human = intervElt.getHuman().get();
        map<string,size_t> identifierMap;
        
        // 1. Read effects
        for( scnXml::HumanInterventions::EffectConstIterator it =
                human.getEffect().begin(), end = human.getEffect().end();
                it != end; ++it )
        {
            const scnXml::HumanInterventionEffect& effect = *it;
            if( identifierMap.count( effect.getId() ) > 0 ){
                ostringstream msg;
                msg << "The id attribute of intervention.human.effect elements must be unique; found \""
                        << effect.getId() << "\" twice.";
                throw util::xml_scenario_error( msg.str() );
            }
            size_t index = humanEffects.size();        // i.e. index of next item
            identifierMap[effect.getId()] = index;
            if( effect.getMDA().present() ){
                //TODO: separate descriptions?
                //TODO: allow continuous deployment
                humanEffects.push_back( new MDAEffect( index, effect.getMDA().get() ) );
            }else if( effect.getPEV().present() ){
                //TODO: allow multiple descriptions of each vaccine type
                humanEffects.push_back( new VaccineEffect( index, effect.getPEV().get(), Host::Vaccine::PEV ) );
            }else if( effect.getBSV().present() ){
                humanEffects.push_back( new VaccineEffect( index, effect.getBSV().get(), Host::Vaccine::BSV ) );
            }else if( effect.getTBV().present() ){
                humanEffects.push_back( new VaccineEffect( index, effect.getTBV().get(), Host::Vaccine::TBV ) );
            }else if( effect.getIPT().present() ){
                humanEffects.push_back( new IPTEffect( index, effect.getIPT().get() ) );
            }else if( effect.getITN().present() ){
                humanEffects.push_back( new ITNEffect( index, effect.getITN().get(), transmission ) );
            }else if( effect.getIRS().present() ){
                humanEffects.push_back( new IRSEffect( index, effect.getIRS().get(), transmission ) );
            }else if( effect.getGVI().present() ){
                if( species_index_map == 0 )
                    species_index_map = &transmission.getSpeciesIndexMap();
                humanEffects.push_back( new interventions::GVIParams( index, effect.getGVI().get(), *species_index_map ) );
            }else if( effect.getCohort().present() ){
                humanEffects.push_back( new CohortSelectionEffect( index ) );
            }else if( effect.getClearImmunity().present() ){
                humanEffects.push_back( new ClearImmunityEffect( index ) );
            }else{
                throw util::xml_scenario_error(
                    "expected intervention.human.effect element to have a "
                    "child, didn't find it (perhaps I need updating)" );
            }
        }
        
        // 2. Read list of interventions
        for( scnXml::HumanInterventions::InterventionConstIterator it =
                human.getIntervention().begin(),
                end = human.getIntervention().end(); it != end; ++it )
        {
            list<Host::Vaccine::Types> vaccineEffects;      // for vaccine EPI deployment
            const scnXml::Intervention& elt = *it;
            // 2.a intervention effects
            HumanIntervention *intervention = new HumanIntervention();
            for( scnXml::Intervention::EffectConstIterator it2 = elt.getEffect().begin(),
                    end2 = elt.getEffect().end(); it2 != end2; ++it2 )
            {
                map<string,size_t>::const_iterator result = identifierMap.find( it2->getId() );
                if( result == identifierMap.end() ){
                    ostringstream msg;
                    msg << "human intervention references effect with id \""
                        << it2->getId()
                        << "\", but no effect with this id was found";
                    throw util::xml_scenario_error( msg.str() );
                }
                const HumanInterventionEffect* effect = &humanEffects[result->second];
                intervention->addEffect( effect );
                const VaccineEffect *vaccEffect = dynamic_cast<const VaccineEffect*>( effect );
                if( vaccEffect != 0 ) {
                    vaccineEffects.push_back( vaccEffect->getVaccineType() );
                }
            }
            intervention->sortEffects();
            
            // 2.b intervention deployments
            if( elt.getContinuous().present() ){
                const scnXml::ContinuousList::DeploySequence& ctsSeq =
                        elt.getContinuous().get().getDeploy();
                for( scnXml::ContinuousList::DeployConstIterator it2 = ctsSeq.begin(),
                    end2 = ctsSeq.end(); it2 != end2; ++it2 )
                {
                    continuous.push_back( new ContinuousHumanIntervention( *it2, intervention ) );
                }
                for( list<Host::Vaccine::Types>::const_iterator it = vaccineEffects.begin(); it != vaccineEffects.end(); ++it ){
                    Host::Vaccine::types[*it].initSchedule( ctsSeq );
                }
            }
            for( scnXml::Intervention::TimedConstIterator timedIt = elt.getTimed().begin();
                timedIt != elt.getTimed().end(); ++timedIt )
            {
                if( timedIt->getCumulativeCoverage().present() ){
                    const scnXml::CumulativeCoverage& cumCov = timedIt->getCumulativeCoverage().get();
                    map<string,size_t>::const_iterator effect_it = identifierMap.find( cumCov.getEffect() );
                    if( effect_it == identifierMap.end() ){
                        ostringstream msg;
                        msg << "interventions.human.intervention.timed."
                                "cumulativeCoverage: element refers to effect identifier \""
                                << cumCov.getEffect()
                                << "\" but no effect with this id was found";
                        throw util::xml_scenario_error( msg.str() );
                    }
                    size_t effect = effect_it->second;
                    TimeStep maxAge = TimeStep::fromYears( cumCov.getMaxAgeYears() );
                    for( scnXml::MassListWithCum::DeployConstIterator it2 =
                            timedIt->getDeploy().begin(), end2 =
                            timedIt->getDeploy().end(); it2 != end2; ++it2 )
                    {
                        timed.push_back( new TimedCumulativeHumanDeployment( *it2, intervention, effect, maxAge ) );
                    }
                }else{
                    for( scnXml::MassListWithCum::DeployConstIterator it2 =
                            timedIt->getDeploy().begin(), end2 =
                            timedIt->getDeploy().end(); it2 != end2; ++it2 )
                    {
                        timed.push_back( new TimedHumanDeployment( *it2, intervention ) );
                    }
                }
            }
            humanInterventions.push_back( intervention );
        }
    }
    if( intervElt.getImportedInfections().present() ){
        const scnXml::ImportedInfections& ii = intervElt.getImportedInfections().get();
        importedInfections.init( ii );
    }
    // Must come after vaccines are initialised:
    if( intervElt.getInsertR_0Case().present() ){
        const scnXml::InsertR_0Case& elt = intervElt.getInsertR_0Case().get();
        if( elt.getTimedDeployment().size() > 0 ){
            Host::Vaccine::verifyEnabledForR_0();
            // timed deployments:
            typedef scnXml::InsertR_0Case::TimedDeploymentSequence::const_iterator It;
            for( It it = elt.getTimedDeployment().begin(); it != elt.getTimedDeployment().end(); ++it ){
                timed.push_back( new TimedR_0Deployment( TimeStep( it->getTime() ) ) );
            }
        }
    }
    if( intervElt.getUninfectVectors().present() ){
        const scnXml::UninfectVectors& elt = intervElt.getUninfectVectors().get();
        if( elt.getTimedDeployment().size() > 0 ){
            // timed deployments:
            typedef scnXml::UninfectVectors::TimedDeploymentSequence::const_iterator It;
            for( It it = elt.getTimedDeployment().begin(); it != elt.getTimedDeployment().end(); ++it ){
                timed.push_back( new TimedUninfectVectorsDeployment( TimeStep( it->getTime() ) ) );
            }
        }
    }
    if( intervElt.getVectorPop().present() ){
        typedef scnXml::VectorPop::InterventionSequence SeqT;
        const SeqT& seq = intervElt.getVectorPop().get().getIntervention();
        size_t instance = 0;
        for( SeqT::const_iterator it = seq.begin(), end = seq.end(); it != end; ++it ){
            const scnXml::VectorIntervention& elt = *it;
            if (elt.getTimed().present() ) {
                population._transmissionModel->initVectorInterv( elt.getDescription().getAnopheles(), instance );
                
                const scnXml::TimedBaseList::DeploySequence& seq = elt.getTimed().get().getDeploy();
                typedef scnXml::TimedBaseList::DeploySequence::const_iterator It;
                for ( It it = seq.begin(); it != seq.end(); ++it ) {
                    timed.push_back( new TimedVectorDeployment(TimeStep( it->getTime() ), instance) );
                }
                instance++;
            }
        }
    }

    // lists must be sorted, increasing
    // For reproducability, we need to use stable_sort, not sort.
    // NOTE: I'd rather use stable_sort, but it's not available. Results are
    // the same without as with a hacked BOOST version including stable_sort.
    continuous.sort();
    timed.sort();
    
    // make sure the list ends with something always in the future, so we don't
    // have to check nextTimed is within range:
    timed.push_back( new DummyTimedDeployment() );
}

void InterventionManager::loadFromCheckpoint( OM::Population& population, TimeStep interventionTime ){
    // We need to re-deploy changeHS and changeEIR interventions, but nothing
    // else. nextTimed should be zero so we can go through all past interventions.
    // Only redeploy those which happened before this timestep.
    assert( nextTimed == 0 );
    while( timed[nextTimed].time < interventionTime ){
        if( dynamic_cast<TimedChangeHSDeployment*>(&timed[nextTimed])!=0 ||
            dynamic_cast<TimedChangeEIRDeployment*>(&timed[nextTimed])!=0 ){
            //Note: neither changeHS nor changeEIR interventions care what the
            //current timestep is when they are deployed, so we don't need to
            //tell them the deployment time.
            timed[nextTimed].deploy( population );
        }
        nextTimed += 1;
    }
}


void InterventionManager::deploy(OM::Population& population) {
    if( TimeStep::interventionPeriod < TimeStep(0) )
        return;
    
    // deploy imported infections (not strictly speaking an intervention)
    importedInfections.import( population );
    
    // deploy timed interventions
    while( timed[nextTimed].time <= TimeStep::interventionPeriod ){
        timed[nextTimed].deploy( population );
        nextTimed += 1;
    }
    
    // deploy continuous interventions
    for( Population::Iter it = population.begin(); it != population.end(); ++it ){
        uint32_t nextCtsDist = it->getNextCtsDist();
        // deploy continuous interventions
        while( nextCtsDist < continuous.size() )
        {
            if( !continuous[nextCtsDist].filterAndDeploy( *it, population ) )
                break;  // deployment (and all remaining) happens in the future
            nextCtsDist = it->incrNextCtsDist();
        }
    }
}

} }
