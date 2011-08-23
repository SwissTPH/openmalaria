/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Transmission/Vector/MosquitoLifeCycle.h"
#include "schema/entomology.h"
#include <cmath>

namespace OM {
namespace Transmission {
namespace Vector {

void MosqLifeCycleParams::initMosqLifeCycle( const scnXml::LifeCycle& lifeCycle ){
    // Simple constants stored in XML:
    eggStageDuration = lifeCycle.getEggStage().getDuration();
    larvalStageDuration = lifeCycle.getLarvalStage().getDuration();
    pupalStageDuration = lifeCycle.getPupalStage().getDuration();
    // we're only interested in female eggs, hence divide by 2:
    fEggsLaidByOviposit = lifeCycle.getEggsLaidByOviposit().getValue() / 2.0;
    //NOTE: store daily or whole-stage probability of survival?
    pSurvEggStage = lifeCycle.getEggStage().getSurvival();
    pSurvDayAsLarvae = std::pow( lifeCycle.getLarvalStage().getSurvival(), 1.0 / larvalStageDuration );
    pSurvPupalStage = lifeCycle.getPupalStage().getSurvival();
    
    // constants varying by larval age; probably stored directly in XML:
    larvaeResourceUsage.reserve( larvalStageDuration );
    effectCompetitionOnLarvae.reserve( larvalStageDuration );
    const scnXml::LarvalStage::DailySequence& larvDev = lifeCycle.getLarvalStage().getDaily();
    for( scnXml::LarvalStage::DailyConstIterator it = larvDev.begin(); it!=larvDev.end(); ++it ){
        larvaeResourceUsage.push_back( it->getResourceUsage() );
        effectCompetitionOnLarvae.push_back( it->getEffectCompetition() );
    }
    
    // complex derivation: annual resource availability to larvae
    //NOTE: this is set by fitLarvalResourcesFromEmergence()
    larvalResources.resize(365);
}


void MosquitoLifeCycle::init( const MosqLifeCycleParams& lcParams ){
    // Shouldn't matter that values start at 0, since the outputs of this model
    // aren't used before all zeros have been replaced.
    newEggs.assign( lcParams.eggStageDuration, 0.0 );
    numLarvae.assign( lcParams.larvalStageDuration, 0.0 );
    newPupae.assign( lcParams.pupalStageDuration, 0.0 );
}

double MosquitoLifeCycle::getResRequirements( const MosqLifeCycleParams& lcParams ) const{
    double resReq = 0.0;
    for( int age=0; age<lcParams.larvalStageDuration; ++age){
        resReq += lcParams.larvaeResourceUsage[age] * numLarvae[age];
    }
    return resReq;
}

double MosquitoLifeCycle::updateEmergence( const MosqLifeCycleParams& lcParams,
                                           double nOvipositingMosqs,
                                           size_t d, size_t dYear1 ){
    // num newly emerging adults comes from num new pupae
    // pupalStageDuration days ago:
    double newAdults =
        lcParams.pSurvPupalStage * newPupae[d % lcParams.pupalStageDuration];
    
    // resource competition during last time-step (L(t) * gamma(t))
    double resourceCompetition = getResRequirements( lcParams )
        * lcParams.larvalResources[dYear1];
    // num new pupae uses larval development formula based on num larvae
    // which were one day away from becoming adults yesterday
    newPupae[d % lcParams.pupalStageDuration] =
        lcParams.pSurvDayAsLarvae * numLarvae[lcParams.larvalStageDuration-1] /
        ( 1.0 + resourceCompetition * lcParams.effectCompetitionOnLarvae[lcParams.larvalStageDuration-1] );
    for( size_t age=lcParams.larvalStageDuration-1;age>=1;--age ){
        numLarvae[age] = lcParams.pSurvDayAsLarvae * numLarvae[age-1] /
        ( 1.0 + resourceCompetition * lcParams.effectCompetitionOnLarvae[age-1] );
    }
    
    // num new larvae comes from num eggs laid eggStageDuration days ago:
    numLarvae[ 0 ] =
        lcParams.pSurvEggStage * newEggs[d % lcParams.eggStageDuration];
    
    // num eggs laid depends on number of mosquitoes which completed a
    // feeding & egg-laying cycle starting tau days ago:
    newEggs[d % lcParams.eggStageDuration] =
        lcParams.fEggsLaidByOviposit * nOvipositingMosqs;
    
    return newAdults;
}

}
}
}
