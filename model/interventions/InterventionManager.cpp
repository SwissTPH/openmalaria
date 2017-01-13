/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

#include "interventions/InterventionManager.hpp"
#include "Population.h"
#include "util/random.h"
#include "util/CommandLine.h"
#include "util/timeConversions.h"
#include "interventions/GVI.h"
#include "interventions/IRS.h"
#include "interventions/ITN.h"
#include "interventions/Vaccine.h"
#include "interventions/HumanInterventionComponents.hpp"
#include "interventions/Deployments.hpp"
#include "WithinHost/Diagnostic.h"

namespace OM { namespace interventions {

// ———  InterventionManager  ———

// static memory:

std::map<std::string,ComponentId> InterventionManager::identifierMap;
boost::ptr_vector<HumanInterventionComponent> InterventionManager::humanComponents;
boost::ptr_vector<HumanIntervention> InterventionManager::humanInterventions;
ptr_vector<ContinuousHumanDeployment> InterventionManager::continuous;
ptr_vector<TimedDeployment> InterventionManager::timed;
uint32_t InterventionManager::nextTimed;
OM::Host::ImportedInfections InterventionManager::importedInfections;

// declared in HumanComponents.h:
vector<ComponentId> removeAtIds[SubPopRemove::NUM];

// static functions:

void InterventionManager::init (const scnXml::Interventions& intervElt){
    nextTimed = 0;
    
    if( intervElt.getChangeHS().present() ){
        const scnXml::ChangeHS& chs = intervElt.getChangeHS().get();
        if( chs.getTimedDeployment().size() > 0 ){
            // timed deployments:
            typedef scnXml::ChangeHS::TimedDeploymentSequence::const_iterator It;
            for( It it = chs.getTimedDeployment().begin(); it != chs.getTimedDeployment().end(); ++it ){
                try{
                    SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                    timed.push_back( new TimedChangeHSDeployment( date, *it ) );
                }catch( const util::format_error& e ){
                    throw util::xml_scenario_error( string("interventions/changeHS/timedDeployment/time: ").append(e.message()) );
                }
            }
        }
    }
    if( intervElt.getChangeEIR().present() ){
        const scnXml::ChangeEIR& eir = intervElt.getChangeEIR().get();
        if( eir.getTimedDeployment().size() > 0 ){
            // timed deployments:
            typedef scnXml::ChangeEIR::TimedDeploymentSequence::const_iterator It;
            for( It it = eir.getTimedDeployment().begin(); it != eir.getTimedDeployment().end(); ++it ){
                try{
                    SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                    timed.push_back( new TimedChangeEIRDeployment( date, *it ) );
                }catch( const util::format_error& e ){
                    throw util::xml_scenario_error( string("interventions/changeEIR/timedDeployment/time: ").append(e.message()) );
                }
            }
        }
    }
    // species_index_map is not available with the non-vector model or
    // non-dynamic mode, so setting it (lazily) also checks sim mode:
    const map<string,size_t>* species_index_map = 0;
    if( intervElt.getHuman().present() ){
        const scnXml::HumanInterventions& human = intervElt.getHuman().get();
        
        // 1. Read components
        for( scnXml::HumanInterventions::ComponentConstIterator it =
                human.getComponent().begin(), end = human.getComponent().end();
                it != end; ++it )
        {
            const scnXml::HumanInterventionComponent& component = *it;
            if( identifierMap.count( component.getId() ) > 0 ){
                ostringstream msg;
                msg << "The id attribute of interventions.human.component elements must be unique; found \""
                        << component.getId() << "\" twice.";
                throw util::xml_scenario_error( msg.str() );
            }
            ComponentId id( humanComponents.size() );        // i.e. index of next item
            identifierMap.insert( make_pair(component.getId(), id) );
            
            SimTime expireAfter = SimTime::future();
            if( component.getSubPopRemoval().present() ){
                const scnXml::SubPopRemoval& removeOpts = component.getSubPopRemoval().get();
                if( removeOpts.getOnFirstBout() ){
                    removeAtIds[SubPopRemove::ON_FIRST_BOUT].push_back( id );
                }
                if( removeOpts.getOnFirstInfection() ){
                    removeAtIds[SubPopRemove::ON_FIRST_INFECTION].push_back( id );
                }
                if( removeOpts.getOnFirstTreatment() ){
                    removeAtIds[SubPopRemove::ON_FIRST_TREATMENT].push_back( id );
                }
                if( removeOpts.getAfterYears().present() ){
                    expireAfter = SimTime::fromYearsN( removeOpts.getAfterYears().get() );
                }
            }
            
            HumanInterventionComponent *hiComponent;
            if( component.getScreen().present() ){
                hiComponent = new ScreenComponent( id, component.getScreen().get() );
            }else if( component.getTreatSimple().present() ){
                hiComponent = new TreatSimpleComponent( id, component.getTreatSimple().get() );
            }else if( component.getTreatPKPD().present() ){
                hiComponent = new TreatPKPDComponent( id, component.getTreatPKPD().get() );
            }else if( component.getDecisionTree().present() ){
                hiComponent = new DecisionTreeComponent( id, component.getDecisionTree().get() );
            }else if( component.getPEV().present() ){
                hiComponent = new VaccineComponent( id, component.getPEV().get(), Vaccine::PEV );
            }else if( component.getBSV().present() ){
                hiComponent = new VaccineComponent( id, component.getBSV().get(), Vaccine::BSV );
            }else if( component.getTBV().present() ){
                hiComponent = new VaccineComponent( id, component.getTBV().get(), Vaccine::TBV );
            }else if( component.getITN().present() ){
                if( species_index_map == 0 )
                    species_index_map = &sim::transmission().getSpeciesIndexMap();
                hiComponent = new ITNComponent( id, component.getITN().get(), *species_index_map );
            }else if( component.getIRS().present() ){
                if( species_index_map == 0 )
                    species_index_map = &sim::transmission().getSpeciesIndexMap();
                hiComponent = new IRSComponent( id, component.getIRS().get(), *species_index_map );
            }else if( component.getGVI().present() ){
                if( species_index_map == 0 )
                    species_index_map = &sim::transmission().getSpeciesIndexMap();
                hiComponent = new GVIComponent( id, component.getGVI().get(), *species_index_map );
            }else if( component.getRecruitmentOnly().present() ){
                hiComponent = new RecruitmentOnlyComponent( id );
            }else if( component.getClearImmunity().present() ){
                hiComponent = new ClearImmunityComponent( id );
            }else{
                throw util::xml_scenario_error(
                    "expected intervention.human.component element to have a "
                    "child, didn't find it (perhaps I need updating)" );
            }
            hiComponent->setExpireAfter( expireAfter );
            humanComponents.push_back( hiComponent );
        }
        
        // 2. Read the list of deployments
        for( scnXml::HumanInterventions::DeploymentConstIterator it =
                human.getDeployment().begin(),
                end = human.getDeployment().end(); it != end; ++it )
        {
            const scnXml::Deployment& elt = *it;
            // 2.a intervention components
            HumanIntervention *intervention = new HumanIntervention( elt.getComponent(),
                elt.getCondition() );
            
            // 2.b intervention deployments
            for( scnXml::Deployment::ContinuousConstIterator ctsIt = elt.getContinuous().begin();
                ctsIt != elt.getContinuous().end(); ++ctsIt )
            {
                ComponentId subPop = ComponentId_pop;
                bool complement = false;
                if( ctsIt->getRestrictToSubPop().present() ){
                    const string& subPopStr = ctsIt->getRestrictToSubPop().get().getId();
                    subPop = getComponentId( subPopStr );
                    complement = ctsIt->getRestrictToSubPop().get().getComplement();
                }
                const scnXml::ContinuousList::DeploySequence& ctsSeq = ctsIt->getDeploy();
                for( scnXml::ContinuousList::DeployConstIterator it2 = ctsSeq.begin(),
                    end2 = ctsSeq.end(); it2 != end2; ++it2 )
                {
                    try{
                        SimTime begin = SimTime::zero();    // intervention period starts at 0
                        if( it2->getBegin().present() ){
                            begin = UnitParse::readDate(it2->getBegin().get(),
                                                        UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        }
                        SimTime end = SimTime::future();
                        if( it2->getEnd().present() ){
                            end = UnitParse::readDate(it2->getEnd().get(),
                                                      UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        }
                        continuous.push_back( new ContinuousHumanDeployment(
                            begin, end, *it2, intervention, subPop, complement ) );
                    }catch( const util::format_error& e ){
                        throw util::xml_scenario_error(
                            string("interventions/human/deployment/continuous/deploy: ")
                            .append(e.message()) );
                    }
                }
            }
            for( scnXml::Deployment::TimedConstIterator timedIt = elt.getTimed().begin();
                timedIt != elt.getTimed().end(); ++timedIt )
            {
                ComponentId subPop = ComponentId_pop;
                bool complement = false;
                if( timedIt->getRestrictToSubPop().present() ){
                    const string& subPopStr = timedIt->getRestrictToSubPop().get().getId();
                    subPop = getComponentId( subPopStr );
                    complement = timedIt->getRestrictToSubPop().get().getComplement();
                }
                try{
                    multimap<SimTime, const scnXml::MassDeployment*> deployTimes;
                    for( scnXml::MassListWithCum::DeployConstIterator it2 =
                            timedIt->getDeploy().begin(), end2 =
                            timedIt->getDeploy().end(); it2 != end2; ++it2 )
                    {
                        SimTime date = UnitParse::readDate(it2->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                        
                        if( it2->getRepeatStep().present() != it2->getRepeatEnd().present() ){
                            throw util::xml_scenario_error( "deploy: use of repeatStep or repeatEnd without other" );
                        }
                        if( it2->getRepeatStep().present() ){
                            SimTime step = UnitParse::readDuration( it2->getRepeatStep().get(), UnitParse::NONE );
                            if( step < SimTime::oneTS() ){
                                throw util::xml_scenario_error( "deploy: repeatStep must be >= 1" );
                            }
                            SimTime end = UnitParse::readDate( it2->getRepeatEnd().get(), UnitParse::NONE );
                            while(date < end){
                                deployTimes.insert(make_pair(date, &*it2));
                                date += step;
                            }
                        }else{
                            deployTimes.insert(make_pair(date, &*it2));
                        }
                    }
                    SimTime lastTime = SimTime::never();
                    for( multimap<SimTime, const scnXml::MassDeployment*>::const_iterator deploy =
                        deployTimes.begin(), dpEnd = deployTimes.end(); deploy != dpEnd; ++deploy )
                    {
                        if( deploy->first == lastTime ){
                            ostringstream msg;
                            msg << "Timed deployment of components ";
                            bool first = true;
                            for( xsd::cxx::tree::sequence<scnXml::Component>::const_iterator cp =
                                elt.getComponent().begin(), cpEnd = elt.getComponent().end(); cp != cpEnd; ++cp )
                            {
                                if( !first ){ msg << ", "; }else{ first=false; }
                                msg << cp->getId();
                            }
                            msg << " has multiple deplyoments at time " << lastTime.date();
                            throw util::xml_scenario_error(msg.str());
                        }
                        lastTime = deploy->first;
                    }
                    
                    if( timedIt->getCumulativeCoverage().present() ){
                        const scnXml::CumulativeCoverage& cumCov = timedIt->getCumulativeCoverage().get();
                        ComponentId cumCovComponent = getComponentId( cumCov.getComponent() );
                        for( multimap<SimTime, const scnXml::MassDeployment*>::const_iterator deploy =
                            deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy )
                        {
                            timed.push_back( new TimedCumulativeHumanDeployment( deploy->first, *deploy->second, intervention, subPop, complement, cumCovComponent ) );
                        }
                    }else{
                        for( multimap<SimTime, const scnXml::MassDeployment*>::const_iterator deploy =
                            deployTimes.begin(), end = deployTimes.end(); deploy != end; ++deploy )
                        {
                            timed.push_back( new TimedHumanDeployment( deploy->first, *deploy->second, intervention, subPop, complement ) );
                        }
                    }
                }catch( const util::format_error& e ){
                    throw util::xml_scenario_error( string("interventions/human/deployment/timed/deploy/time: ").append(e.message()) );
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
        throw util::xml_scenario_error("R_0 code is disabled to reduce "
            "maintenance. If you need it, either use an old OpenMalaria "
            "version (pre 32) or request it be reenabled.");
        // code disabled; search for R_0 to find more blocks like this:
#if 0
        const scnXml::InsertR_0Case& elt = intervElt.getInsertR_0Case().get();
        if( elt.getTimedDeployment().size() > 0 ){
            Vaccine::verifyEnabledForR_0();
            // timed deployments:
            typedef scnXml::InsertR_0Case::TimedDeploymentSequence::const_iterator It;
            for( It it = elt.getTimedDeployment().begin(); it != elt.getTimedDeployment().end(); ++it ){
                SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                timed.push_back( new TimedR_0Deployment( date ) );
            }
        }
#endif
    }
    if( intervElt.getUninfectVectors().present() ){
        const scnXml::UninfectVectors& elt = intervElt.getUninfectVectors().get();
        if( elt.getTimedDeployment().size() > 0 ){
            // timed deployments:
            typedef scnXml::UninfectVectors::TimedDeploymentSequence::const_iterator It;
            for( It it = elt.getTimedDeployment().begin(); it != elt.getTimedDeployment().end(); ++it ){
                SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                timed.push_back( new TimedUninfectVectorsDeployment( date ) );
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
                sim::transmission().initVectorInterv( elt.getDescription().getAnopheles(), instance, elt.getName() );
                
                const scnXml::TimedBaseList::DeploySequence& seq = elt.getTimed().get().getDeploy();
                typedef scnXml::TimedBaseList::DeploySequence::const_iterator It;
                for( It it = seq.begin(); it != seq.end(); ++it ) {
                    SimTime date = UnitParse::readDate(it->getTime(), UnitParse::STEPS /*STEPS is only for backwards compatibility*/);
                    timed.push_back( new TimedVectorDeployment( date, instance ) );
                }
                instance++;
            }
        }
    }
    if( intervElt.getVectorTrap().present() ){
        size_t instance = 0;
        foreach( const scnXml::VectorTrap& trap, intervElt.getVectorTrap().get().getIntervention() ){
            sim::transmission().initVectorTrap(
                    trap.getDescription(), instance, trap.getName() );
            if( trap.getTimed().present() ) {
                foreach( const scnXml::Deploy1 deploy, trap.getTimed().get().getDeploy() ){
                    SimTime time = UnitParse::readDate(deploy.getTime(), UnitParse::STEPS);
                    double ratio = deploy.getRatioToHumans();
                    SimTime lifespan = UnitParse::readDuration(deploy.getLifespan(), UnitParse::NONE);
                    timed.push_back( new TimedTrapDeployment( time, instance, ratio, lifespan ) );
                }
            }
            instance += 1;
        }
    }

    // lists must be sorted, increasing
    // For reproducability, we need to use stable_sort, not sort.
    // NOTE: I'd rather use stable_sort, but it's not available. We also can't
    // use std::stable_sort (from <algorithm>). Results are
    // the same without as with a hacked BOOST version including stable_sort.
    continuous.sort();
    timed.sort();
    
    // make sure the list ends with something always in the future, so we don't
    // have to check nextTimed is within range:
    timed.push_back( new DummyTimedDeployment() );
    
#ifdef WITHOUT_BOINC
    if( util::CommandLine::option( util::CommandLine::PRINT_INTERVENTIONS ) ){
        cout << "Continuous deployments:" << endl
            << "begin\tend\tage\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for( ptr_vector<ContinuousHumanDeployment>::const_iterator it =
            continuous.begin(); it != continuous.end(); ++it ){
            it->print_details( std::cout );
            cout << endl;
        }
        cout << "Timed deployments:" << endl
            << "time\tmin age\tmax age\tsub pop\tcompl\tcoverag\tcomponents" << endl;
        for( ptr_vector<TimedDeployment>::const_iterator it =
            timed.begin(); it != timed.end(); ++it ){
            it->print_details( std::cout );
            cout << endl;
        }
        cout << "Human components:" << endl;
        for( ptr_vector<HumanInterventionComponent>::const_iterator it =
            humanComponents.begin(); it != humanComponents.end(); ++it ){
            it->print_details( cout );
            cout << endl;
        }
    }
#endif
}

ComponentId InterventionManager::getComponentId( const string textId )
{
    map<string,ComponentId>::const_iterator it = identifierMap.find( textId );
    if( it == identifierMap.end() ){
        ostringstream msg;
        msg << "unable to find an intervention component with id \""
            << textId << "\" (wrong name, no definition or used before definition?)";
        throw util::xml_scenario_error( msg.str() );
    }
    return it->second;
}

void InterventionManager::loadFromCheckpoint( SimTime interventionTime ){
    // We need to re-deploy changeHS and changeEIR interventions, but nothing
    // else. nextTimed should be zero so we can go through all past interventions.
    // Only redeploy those which happened before this time step.
    assert( nextTimed == 0 );
    while( timed[nextTimed].time < interventionTime ){
        TimedDeployment *deployment = &timed[nextTimed];
        if( dynamic_cast<TimedChangeHSDeployment*>(deployment)!=0 ||
            dynamic_cast<TimedChangeEIRDeployment*>(deployment)!=0 ){
            //Note: neither changeHS nor changeEIR interventions care what the
            //current time step is when they are deployed, so we don't need to
            //tell them the deployment time.
            deployment->deploy( sim::humanPop() );
        }
        nextTimed += 1;
    }
}


void InterventionManager::deploy(OM::Population& population) {
    if( sim::intervNow() < SimTime::zero() )
        return;
    
    // deploy imported infections (not strictly speaking an intervention)
    importedInfections.import( population );
    
    // deploy timed interventions
    while( timed[nextTimed].time <= sim::intervNow() ){
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
