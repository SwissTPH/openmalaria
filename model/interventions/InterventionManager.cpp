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

#include "interventions/InterventionManager.hpp"
#include "Population.h"
#include "util/random.h"
#include "util/CommandLine.h"
#include "Monitoring/Survey.h"
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
    
boost::ptr_vector<HumanInterventionComponent> InterventionManager::humanComponents;
boost::ptr_vector<HumanIntervention> InterventionManager::humanInterventions;
ptr_vector<ContinuousHumanDeployment> InterventionManager::continuous;
ptr_vector<TimedDeployment> InterventionManager::timed;
uint32_t InterventionManager::nextTimed;
OM::Host::ImportedInfections InterventionManager::importedInfections;
bool InterventionManager::_cohortEnabled;

// static functions:

void InterventionManager::init (const scnXml::Interventions& intervElt, OM::Population& population){
    nextTimed = 0;
    _cohortEnabled = false;
    
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
        map<string,ComponentId> identifierMap;
        
        // 1. Read components
        for( scnXml::HumanInterventions::ComponentConstIterator it =
                human.getComponent().begin(), end = human.getComponent().end();
                it != end; ++it )
        {
            const scnXml::HumanInterventionComponent& component = *it;
            if( identifierMap.count( component.getId() ) > 0 ){
                ostringstream msg;
                msg << "The id attribute of intervention.human.component elements must be unique; found \""
                        << component.getId() << "\" twice.";
                throw util::xml_scenario_error( msg.str() );
            }
            ComponentId id( humanComponents.size() );        // i.e. index of next item
            identifierMap.insert( make_pair(component.getId(), id) );
            if( component.getMDA().present() ){
                //TODO(monitoring): report
                humanComponents.push_back( new MDAComponent( id, component.getMDA().get() ) );
            }else if( component.getMDA1D().present() ){
                //TODO(monitoring): report
                humanComponents.push_back( new MDA1DComponent( id, component.getMDA1D().get() ) );
            }else if( component.getPEV().present() ){
                humanComponents.push_back( new VaccineComponent( id, component.getPEV().get(), Vaccine::PEV ) );
            }else if( component.getBSV().present() ){
                humanComponents.push_back( new VaccineComponent( id, component.getBSV().get(), Vaccine::BSV ) );
            }else if( component.getTBV().present() ){
                humanComponents.push_back( new VaccineComponent( id, component.getTBV().get(), Vaccine::TBV ) );
            }else if( component.getITN().present() ){
                if( species_index_map == 0 )
                    species_index_map = &transmission.getSpeciesIndexMap();
                humanComponents.push_back( new ITNComponent( id, component.getITN().get(), *species_index_map ) );
            }else if( component.getIRS().present() ){
                if( species_index_map == 0 )
                    species_index_map = &transmission.getSpeciesIndexMap();
                humanComponents.push_back( new IRSComponent( id, component.getIRS().get(), *species_index_map ) );
            }else if( component.getGVI().present() ){
                if( species_index_map == 0 )
                    species_index_map = &transmission.getSpeciesIndexMap();
                humanComponents.push_back( new GVIComponent( id, component.getGVI().get(), *species_index_map ) );
            }else if( component.getCohort().present() ){
                humanComponents.push_back( new CohortSelectionComponent( id, component.getCohort().get() ) );
                _cohortEnabled = true;
            }else if( component.getClearImmunity().present() ){
                humanComponents.push_back( new ClearImmunityComponent( id ) );
            }else{
                throw util::xml_scenario_error(
                    "expected intervention.human.component element to have a "
                    "child, didn't find it (perhaps I need updating)" );
            }
        }
        
        // 2. Read the list of deployments
        for( scnXml::HumanInterventions::DeploymentConstIterator it =
                human.getDeployment().begin(),
                end = human.getDeployment().end(); it != end; ++it )
        {
            const scnXml::Deployment& elt = *it;
            // 2.a intervention components
            HumanIntervention *intervention = new HumanIntervention();
            for( scnXml::Deployment::ComponentConstIterator it2 = elt.getComponent().begin(),
                    end2 = elt.getComponent().end(); it2 != end2; ++it2 )
            {
                map<string,ComponentId>::const_iterator result = identifierMap.find( it2->getId() );
                if( result == identifierMap.end() ){
                    ostringstream msg;
                    msg << "human intervention references component with id \""
                        << it2->getId()
                        << "\", but no component with this id was found";
                    throw util::xml_scenario_error( msg.str() );
                }
                const HumanInterventionComponent* component = &humanComponents[result->second.id];
                intervention->addComponent( component );
            }
            intervention->sortComponents();
            
            // 2.b intervention deployments
            for( scnXml::Deployment::ContinuousConstIterator ctsIt = elt.getContinuous().begin();
                ctsIt != elt.getContinuous().end(); ++ctsIt )
            {
                ComponentId cohort = ComponentId_pop;
                if( ctsIt->getRestrictToCohort().present() ){
                    const string& cohortName = ctsIt->getRestrictToCohort().get().getId();
                    map<string,ComponentId>::const_iterator comp_it = identifierMap.find( cohortName );
                    if( comp_it == identifierMap.end() ){
                        ostringstream msg;
                        msg << "interventions.human.intervention.continuous."
                                "restrictToCohort: element refers to cohort "
                                "component with identifier \"" << cohortName
                                << "\" but no component with this id was found";
                        throw util::xml_scenario_error( msg.str() );
                    }
                    cohort = comp_it->second;
                }
                const scnXml::ContinuousList::DeploySequence& ctsSeq = ctsIt->getDeploy();
                for( scnXml::ContinuousList::DeployConstIterator it2 = ctsSeq.begin(),
                    end2 = ctsSeq.end(); it2 != end2; ++it2 )
                {
                    continuous.push_back( new ContinuousHumanDeployment( *it2, intervention, cohort ) );
                }
            }
            for( scnXml::Deployment::TimedConstIterator timedIt = elt.getTimed().begin();
                timedIt != elt.getTimed().end(); ++timedIt )
            {
                ComponentId cohort = ComponentId_pop;
                if( timedIt->getRestrictToCohort().present() ){
                    const string& cohortName = timedIt->getRestrictToCohort().get().getId();
                    map<string,ComponentId>::const_iterator comp_it = identifierMap.find( cohortName );
                    if( comp_it == identifierMap.end() ){
                        ostringstream msg;
                        msg << "interventions.human.intervention.timed."
                                "restrictToCohort: element refers to cohort "
                                "component with identifier \"" << cohortName
                                << "\" but no component with this id was found";
                        throw util::xml_scenario_error( msg.str() );
                    }
                    cohort = comp_it->second;
                }
                if( timedIt->getCumulativeCoverage().present() ){
                    const scnXml::CumulativeCoverage& cumCov = timedIt->getCumulativeCoverage().get();
                    map<string,ComponentId>::const_iterator comp_it = identifierMap.find( cumCov.getComponent() );
                    if( comp_it == identifierMap.end() ){
                        ostringstream msg;
                        msg << "interventions.human.intervention.timed."
                                "cumulativeCoverage: element refers to component identifier \""
                                << cumCov.getComponent()
                                << "\" but no component with this id was found";
                        throw util::xml_scenario_error( msg.str() );
                    }
                    ComponentId component = comp_it->second;
                    TimeStep maxAge = TimeStep::fromYears( cumCov.getMaxAgeYears() );
                    for( scnXml::MassListWithCum::DeployConstIterator it2 =
                            timedIt->getDeploy().begin(), end2 =
                            timedIt->getDeploy().end(); it2 != end2; ++it2 )
                    {
                        timed.push_back( new TimedCumulativeHumanDeployment( *it2, intervention, cohort, component, maxAge ) );
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
                timed.push_back( new TimedR_0Deployment( TimeStep( it->getTime() ) ) );
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
                population._transmissionModel->initVectorInterv( elt.getDescription().getAnopheles(), instance, elt.getName() );
                
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
            << "begin\tend\tage\tcohort\tcoverag\tcomponents" << endl;
        for( ptr_vector<ContinuousHumanDeployment>::const_iterator it =
            continuous.begin(); it != continuous.end(); ++it ){
            it->print_details( std::cout );
            cout << endl;
        }
        cout << "Timed deployments:" << endl
            << "time\tmin age\tmax age\tcohort\tcoverag\tcomponents" << endl;
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
