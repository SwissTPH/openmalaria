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
#include "Population.h"
#include "util/random.h"
#include "Monitoring/Surveys.h"
#include "interventions/GVI.h"
#include "interventions/IRS.h"
#include "interventions/ITN.h"
#include "interventions/HumanInterventionEffects.hpp"
#include "interventions/TimedDeployments.hpp"

namespace OM { namespace interventions {

// ———  ContinuousHumanDeployment  ———

/** Interface for continuous deployment of an intervention. */
class ContinuousHumanDeployment {
public:
    /// Create, passing deployment age
    ContinuousHumanDeployment( const ::scnXml::ContinuousDeployment& elt,
                                 const HumanIntervention* intervention, size_t cohort ) :
            begin( elt.getBegin() ),
            end( elt.getEnd() ),
            deployAge( TimeStep::fromYears( elt.getTargetAgeYrs() ) ),
            cohort( cohort ),
            coverage( elt.getCoverage() ),
            intervention( intervention )
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
    
    /// For sorting
    inline bool operator<( const ContinuousHumanDeployment& that )const{
        return this->deployAge < that.deployAge;
    }
    
    /** Apply filters and potentially deploy.
     * 
     * @returns false iff this deployment (and thus all later ones in the
     *  ordered list) happens in the future. */
    bool filterAndDeploy( Host::Human& human, const Population& population ) const{
        TimeStep age = TimeStep::simulation - human.getDateOfBirth();
        if( deployAge > age ){
            // stop processing continuous deployments for this
            // human for now because remaining ones happen in the future
            return false;
        }else if( deployAge == age ){
            if( begin <= TimeStep::interventionPeriod &&
                TimeStep::interventionPeriod < end &&
                ( human.isInCohort( cohort ) ) &&
                util::random::uniform_01() < coverage )     // RNG call should be last test
            {
                deploy( human, population );
            }
        }//else: for some reason, a deployment age was missed; ignore it
        return true;
    }
    
protected:
    /// Deploy to a selected human.
    void deploy( Host::Human& human, const Population& population ) const{
        intervention->deploy( human, Deployment::CTS );
    }
    
    TimeStep begin, end;    // first timeStep active and first timeStep no-longer active
    TimeStep deployAge;
    size_t cohort;      // size_t maximum value if no cohort
    double coverage;
    const HumanIntervention *intervention;
};


// ———  HumanInterventionEffect  ———

void HumanIntervention::deploy( Human& human, Deployment::Method method ) const{
    for( vector<const HumanInterventionEffect*>::const_iterator it = effects.begin();
            it != effects.end(); ++it )
    {
        const interventions::HumanInterventionEffect& effect = **it;
        effect.deploy( human, method );
        human.updateLastDeployed( effect.getIndex() );
    }
}

bool effectCmp(const HumanInterventionEffect *a, const HumanInterventionEffect *b){
    return a->effectType() < b->effectType();
}
void HumanIntervention::sortEffects(){
    std::stable_sort( effects.begin(), effects.end(), effectCmp );
}


// ———  InterventionManager  ———

auto_ptr<InterventionManager> InterventionManager::instance;

void InterventionManager::makeInstance( const scnXml::Interventions& intervElt, OM::Population& population ){
    instance = auto_ptr<InterventionManager>( new InterventionManager( intervElt, population ) );
}

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
    // species_index_map is not available with the non-vector model or
    // non-dynamic mode, so setting it (lazily) also checks sim mode:
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
                //TODO: monitoring
                humanEffects.push_back( new MDAEffect( index, effect.getMDA().get() ) );
            }else if( effect.getMDA1D().present() ){
                //TODO: monitoring
                humanEffects.push_back( new MDA1DEffect( index, effect.getMDA1D().get() ) );
            }else if( effect.getPEV().present() ){
                //TODO: allow multiple descriptions of each vaccine type
                humanEffects.push_back( new VaccineEffect( index, effect.getPEV().get(), Vaccine::PEV ) );
            }else if( effect.getBSV().present() ){
                humanEffects.push_back( new VaccineEffect( index, effect.getBSV().get(), Vaccine::BSV ) );
            }else if( effect.getTBV().present() ){
                humanEffects.push_back( new VaccineEffect( index, effect.getTBV().get(), Vaccine::TBV ) );
            }else if( effect.getIPT().present() ){
		// TODO: also remove XML elements from XSD in a later versions
                throw util::xml_scenario_error( "The IPT model is no longer available. Use MDA instead." );
            }else if( effect.getITN().present() ){
                if( species_index_map == 0 )
                    species_index_map = &transmission.getSpeciesIndexMap();
                humanEffects.push_back( new ITNEffect( index, effect.getITN().get(), *species_index_map ) );
            }else if( effect.getIRS().present() ){
                if( species_index_map == 0 )
                    species_index_map = &transmission.getSpeciesIndexMap();
                humanEffects.push_back( new IRSEffect( index, effect.getIRS().get(), *species_index_map ) );
            }else if( effect.getGVI().present() ){
                if( species_index_map == 0 )
                    species_index_map = &transmission.getSpeciesIndexMap();
                humanEffects.push_back( new GVIEffect( index, effect.getGVI().get(), *species_index_map ) );
            }else if( effect.getCohort().present() ){
                humanEffects.push_back( new CohortSelectionEffect( index, effect.getCohort().get() ) );
                _cohortEnabled = true;
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
            list<Vaccine::Types> vaccineEffects;      // for vaccine EPI deployment
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
            for( scnXml::Intervention::ContinuousConstIterator ctsIt = elt.getContinuous().begin();
                ctsIt != elt.getContinuous().end(); ++ctsIt )
            {
                size_t cohort = numeric_limits<size_t>::max();
                if( ctsIt->getRestrictToCohort().present() ){
                    const string& id = ctsIt->getRestrictToCohort().get().getId();
                    map<string,size_t>::const_iterator effect_it = identifierMap.find( id );
                    if( effect_it == identifierMap.end() ){
                        ostringstream msg;
                        msg << "interventions.human.intervention.continuous."
                                "restrictToCohort: element refers to cohort "
                                "effect with identifier \"" << id
                                << "\" but no effect with this id was found";
                        throw util::xml_scenario_error( msg.str() );
                    }
                    cohort = effect_it->second;
                }
                const scnXml::ContinuousList::DeploySequence& ctsSeq = ctsIt->getDeploy();
                for( scnXml::ContinuousList::DeployConstIterator it2 = ctsSeq.begin(),
                    end2 = ctsSeq.end(); it2 != end2; ++it2 )
                {
                    continuous.push_back( new ContinuousHumanDeployment( *it2, intervention, cohort ) );
                }
                for( list<Vaccine::Types>::const_iterator it = vaccineEffects.begin(); it != vaccineEffects.end(); ++it ){
                    Vaccine::initSchedule( *it, ctsSeq );
                }
            }
            for( scnXml::Intervention::TimedConstIterator timedIt = elt.getTimed().begin();
                timedIt != elt.getTimed().end(); ++timedIt )
            {
                size_t cohort = numeric_limits<size_t>::max();
                if( timedIt->getRestrictToCohort().present() ){
                    const string& id = timedIt->getRestrictToCohort().get().getId();
                    map<string,size_t>::const_iterator effect_it = identifierMap.find( id );
                    if( effect_it == identifierMap.end() ){
                        ostringstream msg;
                        msg << "interventions.human.intervention.timed."
                                "restrictToCohort: element refers to cohort "
                                "effect with identifier \"" << id
                                << "\" but no effect with this id was found";
                        throw util::xml_scenario_error( msg.str() );
                    }
                    cohort = effect_it->second;
                }
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
                        timed.push_back( new TimedCumulativeHumanDeployment( *it2, intervention, cohort, effect, maxAge ) );
                    }
                }else{
                    for( scnXml::MassListWithCum::DeployConstIterator it2 =
                            timedIt->getDeploy().begin(), end2 =
                            timedIt->getDeploy().end(); it2 != end2; ++it2 )
                    {
                        timed.push_back( new TimedHumanDeployment( *it2, intervention, cohort ) );
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
            Vaccine::verifyEnabledForR_0();
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
