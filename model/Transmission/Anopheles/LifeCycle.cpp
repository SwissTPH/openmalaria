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

#include "Transmission/Anopheles/LifeCycle.h"
#include "schema/entomology.h"
#include <cmath>

namespace OM {
namespace Transmission {
namespace Anopheles {

void LifeCycleParams::initLifeCycle( const scnXml::LifeCycle& lifeCycle ){
    // Simple constants stored in XML:
    eggStageDuration = lifeCycle.getEggStage().getDuration();
    larvalStageDuration = lifeCycle.getLarvalStage().getDuration();
    pupalStageDuration = lifeCycle.getPupalStage().getDuration();
    // we're only interested in female eggs, hence divide by 2:
    fEggsLaidByOviposit = lifeCycle.getFemaleEggsLaidByOviposit().getValue();
    //NOTE: store daily or whole-stage probability of survival?
    pSurvEggStage = lifeCycle.getEggStage().getSurvival();
    pSurvDayAsLarvae = std::pow( lifeCycle.getLarvalStage().getSurvival(), 1.0 / larvalStageDuration );
    pSurvPupalStage = lifeCycle.getPupalStage().getSurvival();
    estimatedLarvalResources = lifeCycle.getEstimatedLarvalResources();
    
    // constants varying by larval age; probably stored directly in XML:
    larvaeResourceUsage.reserve( larvalStageDuration );
    effectCompetitionOnLarvae.reserve( larvalStageDuration );
    const scnXml::LarvalStage::DailySequence& larvDev = lifeCycle.getLarvalStage().getDaily();
    for( scnXml::LarvalStage::DailyConstIterator it = larvDev.begin(); it!=larvDev.end(); ++it ){
        larvaeResourceUsage.push_back( it->getResourceUsage() );
        effectCompetitionOnLarvae.push_back( it->getEffectCompetition() );
    }
    
    // This is set later, by ResourceFitter class.
    invLarvalResources.resize(sim::oneYear());
}

double LifeCycleParams::getResAvailability() const{
    double val = 0.0;
    // get resources available on the next time step
    for( SimTime d = sim::now(), end = sim::now() + sim::oneTS(); d < end;
            d += sim::oneDay() ){
        val += 1.0 / invLarvalResources[mod_nn(d, sim::oneYear())];
    }
    //TODO: why are units per time step?
    return val / sim::oneTS().inDays();
}

void LifeCycle::init( const LifeCycleParams& lcParams ){
    // Shouldn't matter that values start at 0, since the outputs of this model
    // aren't used before all zeros have been replaced.
    newEggs.assign( lcParams.eggStageDuration, 0.0 );
    numLarvae.assign( lcParams.larvalStageDuration, 0.0 );
    newPupae.assign( lcParams.pupalStageDuration, 0.0 );
}

double LifeCycle::getResRequirements( const LifeCycleParams& lcParams ) const{
    double resReq = 0.0;
    for( int age=0; age<lcParams.larvalStageDuration; ++age){
        resReq += lcParams.larvaeResourceUsage[age] * numLarvae[age];
    }
    return resReq;
}

double LifeCycle::updateEmergence( const LifeCycleParams& lcParams,
                                           double nOvipositingMosqs,
                                           SimTime d0 ){
    SimTime d1 = d0 + sim::oneDay();
    // Day of year. Note that d==sim::oneTS() corresponds to Jan 1st, index 0.
    SimTime dYear1 = mod_nn(d0, sim::oneYear());
    // num newly emerging adults comes from num new pupae pupalStageDuration days ago:
    size_t pupaeIndex = mod_nn(d1, lcParams.pupalStageDuration);
    double newAdults = lcParams.pSurvPupalStage * newPupae[pupaeIndex];
    
    // resource competition during last time-step (L(t) * gamma(t))
    double resourceCompetition = getResRequirements( lcParams )
        * lcParams.invLarvalResources[dYear1];  // TODO: scale for larviciding here
    // num new pupae uses larval development formula based on num larvae
    // which were one day away from becoming adults yesterday
    newPupae[pupaeIndex] =
        lcParams.pSurvDayAsLarvae * numLarvae[lcParams.larvalStageDuration-1] /
        ( 1.0 + resourceCompetition * lcParams.effectCompetitionOnLarvae[lcParams.larvalStageDuration-1] );
    for( size_t age=lcParams.larvalStageDuration-1;age>=1;--age ){
        numLarvae[age] = lcParams.pSurvDayAsLarvae * numLarvae[age-1] /
        ( 1.0 + resourceCompetition * lcParams.effectCompetitionOnLarvae[age-1] );
    }
    
    // num new larvae comes from num eggs laid eggStageDuration days ago:
    size_t eggIndex = mod_nn(d1, lcParams.eggStageDuration);
    numLarvae[ 0 ] = lcParams.pSurvEggStage * newEggs[eggIndex];
    
    // num eggs laid depends on number of mosquitoes which completed a
    // feeding & egg-laying cycle starting tau days ago:
    newEggs[eggIndex] = lcParams.fEggsLaidByOviposit * nOvipositingMosqs;
    
    return newAdults;
}

}
}
}
