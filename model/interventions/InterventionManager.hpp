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
#ifndef OM_INTERVENTIONS_INTERVENTIONS
#define OM_INTERVENTIONS_INTERVENTIONS

#include "Global.h"
#include "interventions/Interfaces.hpp"
#include "Host/ImportedInfections.h"
#include "schema/interventions.h"
#include <boost/ptr_container/ptr_vector.hpp>

namespace OM {
    class Population;
    namespace Host {
        class Human;
    }

namespace interventions {
    using ::boost::ptr_vector;

class ContinuousHumanDeployment;
class TimedDeployment;

/** Management of interventions deployed on a per-time-step basis. */
class InterventionManager {
public:
    /** Read XML descriptions. */
    static void init(const scnXml::Interventions& intervElt);
    
    /// Checkpointing
    template<class S>
    static void checkpoint (S& stream) {
        using namespace OM::util::checkpoint;
        // most members are only set from XML,
        // nextTimed varies but is re-set by loadFromCheckpoint
        importedInfections & stream;
    }

    /** Call after loading a checkpoint, passing the intervention-period time.
     * 
     * Serves to replace health-system and EIR where changeHS/changeEIR
     * interventions have been used. */
    static void loadFromCheckpoint( SimTime interventionTime );
    
    /** @brief Deploy interventions
     *
     * Timed interventions are deployed for this time step.
     * 
     * Continuous interventions are deployed as humans reach the target ages.
     * Unlike with vaccines, missing one schedule doesn't preclude the next. */
    static void deploy (OM::Population& population);
    
    /** Get a constant reference to a component class with a certain index.
     * 
     * @throws util::base_exception if the index is out-of-range */
    inline static const HumanInterventionComponent& getComponent( ComponentId id ){
        if( id.id >= humanComponents.size() )
            throw util::base_exception( "invalid component id" );
        return humanComponents[id.id];
    }
    
    /** Get a numeric ComponentId from the textual identifier used in the XML.
     * 
     * If textId is unknown, an xml_scenario_error is thrown. */
    static ComponentId getComponentId( const std::string textId );
    
private:
    // Map of textual identifiers to numeric identifiers for components
    static std::map<std::string,ComponentId> identifierMap;
    // TODO: use C++11 move semantics
    // All human intervention components, indexed by a number. This list is used
    // during initialisation and thereafter only for memory management.
    static boost::ptr_vector<HumanInterventionComponent> humanComponents;
    // All human interventions. These are stored here for memory management
    // only (so that they are deleted when this class is destroyed).
    static boost::ptr_vector<HumanIntervention> humanInterventions;
    // Continuous interventions, sorted by deployment age (weakly increasing)
    static ptr_vector<ContinuousHumanDeployment> continuous;
    // List of all timed interventions. Should be sorted (time weakly increasing).
    static ptr_vector<TimedDeployment> timed;
    static uint32_t nextTimed;  // not chcekpointed (see loadFromCheckpoint)
    
    // imported infections are not really interventions, and handled by a separate class
    // (but are grouped here for convenience and due toassociation in schema)
    static OM::Host::ImportedInfections importedInfections;
};

} }
#endif
