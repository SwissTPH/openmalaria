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

#include "WithinHost/WHVivax.h"
#include "util/random.h"
#include "util/errors.h"
#include "inputData.h"
#include <algorithm>
#include <limits>
#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;

// ———  parameters  ———

// Set from the parameters block:
TimeStep latentp;       // attribute on parameters block

//TODO: these parameters need to be initialised from XML, and possibly made
// static members of VivaxInfection / WHVivax

// current values need replacing:
double probPatentInfectMosquitoes = 1;
double pHetNoPQ = 0;
double pReceivePQ = 0;
double pPqIsEffective = 0;

// Parameters from Amanda:
const int maxNumberHypnozoites = 11;  // chosen for convenience
//TODO: all below need to go in XML
const double baseNumberHypnozoites = 0.8;
const double muReleaseHypnozoite = 2.92;
const double sigmaReleaseHypnozoite = 0.956;
const TimeStep minReleaseHypnozoite( 2 );
const TimeStep bloodStageProtectionLatency( 2 );
double bloodStageLengthWeibullScale = 23.01164;
double bloodStageLengthWeibullShape = 2.331;


// ———  individual models  ———

// number of hypnozoites per brood:
map<double,int> nHypnozoitesProbMap;
void initNHypnozoites(){
    double total = 0.0;
    for( size_t n = 0; n <= maxNumberHypnozoites; ++n )
        total += pow( baseNumberHypnozoites, n );
    
    double cumP = 0.0;
    for( size_t n = 0; n <= maxNumberHypnozoites; ++n ){
        cumP += pow( baseNumberHypnozoites, n ) / total;
        // pair n with the cumulative probability of sampling n:
        nHypnozoitesProbMap[cumP] = n;
    }
}
int sampleNHypnozoites(){
    double x = util::random::uniform_01();
    // upper_bound finds the first key (cumulative probability) greater than x:
    map<double,int>::const_iterator it = nHypnozoitesProbMap.upper_bound( x );
    assert( it != nHypnozoitesProbMap.end() );  // i.e. that we did find a real key
    return it->second;  // corresponding n
}

// time to hypnozoite release after initial release:
TimeStep sampleReleaseDelay(){
    TimeStep delay;
    do{
        delay = TimeStep::fromNearest( util::random::log_normal( muReleaseHypnozoite, sigmaReleaseHypnozoite ) );
    }while( delay < minReleaseHypnozoite );
    return delay;
}


// ———  per-brood code  ———

struct TimeStepCompReverse{
    bool operator()( const TimeStep& a, const TimeStep& b ){
        return b < a;
    }
} timeStepCompR;
VivaxBrood::VivaxBrood(){    
    // primary blood stage plus hypnozoites (relapses)
    releaseDates.push_back( TimeStep::simulation + latentp );
    int numberHypnozoites = sampleNHypnozoites();
    for( int i = 0; i < numberHypnozoites; ++i ){
        TimeStep timeToRelease = TimeStep::simulation + latentp + sampleReleaseDelay();
        releaseDates.push_back( TimeStep::simulation + timeToRelease );
    }
    
    // Sort by time to release, smallest (soonest) last. Explanation of code:
    // http://www.cplusplus.com/reference/algorithm/sort/
    sort( releaseDates.begin(), releaseDates.end(), timeStepCompR );
}

bool VivaxBrood::update( bool& anyNewBloodStage ){
    if( bloodStageClearDate == TimeStep::simulation ){
        //NOTE: this effectively means that both asexual and sexual stage
        // parasites self-terminate. It also means the immune system can
        // protect against new blood-stage infections for a short time.
    }
    
    while( releaseDates.back() == TimeStep::simulation ){
        releaseDates.pop_back();
        // an existing or recently terminated blood stage from the same brood
        // protects against a newly released Hypnozoite
        //NOTE: this is an immunity effect: should there be no immunity when a blood stage first emerges?
        if( bloodStageClearDate + bloodStageProtectionLatency
            >= TimeStep::simulation ) continue;
        
        anyNewBloodStage = true;
        double lengthDays = random::weibull( bloodStageLengthWeibullScale, bloodStageLengthWeibullShape );
        bloodStageClearDate = TimeStep::simulation + TimeStep::fromDaysNearest( lengthDays );
        // Assume gametocytes emerge at the same time (they mature quickly and
        // we have little data, thus assume coincedence of start)
    }
    
    bool isFinished = (!isPatent()) && (releaseDates.size() == 0);
    return isFinished;
}

void VivaxBrood::treatmentBS(){
    // Blood stage treatment: clear both asexual and sexual parasites from the
    // blood. NOTE: we assume infections removed via treatment do not leave
    // protective immunity since the patient was unable to self-clear.
    bloodStageClearDate = TimeStep::never;
}

void VivaxBrood::treatmentLS(){
    releaseDates.clear();       // 100% clearance
    
    /* partial clearance code, in case of need:
    vector<TimeStep> survivingZoites;
    survivingZoites.reserve( releaseDates.size() );   // maximum size we need
    for( vector<TimeStep>::const_iterator it = releaseDates.begin(); it != releaseDates.end(); ++it ){
        if( !random::bernoulli( pClearEachHypnozoite ) ){
            survivingZoites.push_back( *it );    // copy            
        }
    }
    swap( releaseDates, survivingZoites );
    */
}


// ———  per-host code  ———

WHVivax::WHVivax(){
    noPQ = ( pHetNoPQ > 0.0 && random::bernoulli(pHetNoPQ) );
}

void WHVivax::setComorbidityFactor(double factor){
    if( factor != 1.0 )
        throw TRACED_EXCEPTION_DEFAULT( "This vivax model cannot be used with comorbidity heterogeneity" );
}

WHVivax::~WHVivax(){
}

double WHVivax::probTransmissionToMosquito(TimeStep ageTimeSteps, double tbvEfficacy) const{
    for (std::list<VivaxBrood>::const_iterator inf = infections.begin();
         inf != infections.end(); ++inf)
    {
        if( inf->isPatent() ){
            // we have gametocytes from at least one brood
            return probPatentInfectMosquitoes * (1.0 - tbvEfficacy );
        }
    }
    return 0;   // no gametocytes
}

bool WHVivax::summarize(Monitoring::Survey& survey, Monitoring::AgeGroup ageGroup){
    InfectionCount count = countInfections();
    if( count.total != 0 ){
        survey.reportInfectedHosts(ageGroup, 1);
        survey.addToInfections(ageGroup, count.total);
        survey.addToPatentInfections(ageGroup, count.patent);
    }
    if( count.patent > 0 ){
        survey.reportPatentHosts(ageGroup, 1);
        return true;
    }
    return false;
}

void WHVivax::importInfection(){
    // this means one new liver stage infection, which can result in multiple blood stages
    infections.resize( infections.size() + 1 );
}

void WHVivax::update(int nNewInfs, double ageInYears, double BSVEfficacy){
    // create new infections, letting the constructor do the initialisation work:
    infections.resize( infections.size() + nNewInfs );
    
    // update infections
    // NOTE: currently no BSV model
    bool anyNewBloodStage = false;
    std::list<VivaxBrood>::iterator inf = infections.begin();
    while( inf != infections.end() ){
        bool isFinished = inf->update( anyNewBloodStage );
        if( isFinished ) inf = infections.erase( inf );
        else ++inf;
    }
    
    morbidity = Pathogenesis::NONE;
    //NOTE: released hypnozoites blocked by immunity also don't cause fever
    if( anyNewBloodStage ){
        //TODO: 3 probabilities: for initial infection, 0 for relapse given no initial fever,
        // another probability given initial fever
        //TODO: also depends on number of broods ever seen by host (as with cumulativeh for Pf model)
        //NOTE: currently we don't model co-infection or indirect deaths
        /*
        double x = random::uniform_01();
        if( x < pSevMalaria )
            morbidity = Pathogenesis::STATE_SEVERE;
        else if( x < pSevPlusUC )
            morbidity = Pathogenesis::STATE_MALARIA;
        else if( x < pSevPlusUCPlusNMF )
            morbidity = Pathogenesis::STATE_NMF;
        */
    }
}

bool WHVivax::diagnosticDefault() const{
    for (std::list<VivaxBrood>::const_iterator inf = infections.begin();
         inf != infections.end(); ++inf)
    {
        if (inf->isPatent())
            return true;        // at least one patent infection
    }
    return false;
}
bool WHVivax::diagnosticMDA() const{
    return diagnosticDefault();
}

Pathogenesis::State WHVivax::determineMorbidity(double ageYears){
    return morbidity;
}

void WHVivax::immuneSuppression(){
    throw TRACED_EXCEPTION_DEFAULT( "vivax model does not include immune suppression" );
}

WHInterface::InfectionCount WHVivax::countInfections () const{
    InfectionCount count;       // constructor initialises counts to 0
    count.total = infections.size();
    for (std::list<VivaxBrood>::const_iterator inf = infections.begin(); inf != infections.end(); ++inf) {
        if (inf->isPatent())
            count.patent += 1;
    }
    return count;
}
void WHVivax::effectiveTreatment(){
    // This means clear blood stage infection(s) but not liver stage.
    for( list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
        it->treatmentBS();
    }
}
bool WHVivax::optionalPqTreatment(){
    // PQ clears liver stages. We don't worry about the effect of PQ on
    // gametocytes, because these are always cleared by blood-stage drugs with
    // Vivax, and PQ is not given without BS drugs. NOTE: this ignores drug failure.
    if (pReceivePQ > 0.0 && !noPQ && random::bernoulli(pReceivePQ)){
        if( random::bernoulli(pPqIsEffective) ){
            for( list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
                it->treatmentLS();
            }
        }
        return true;    // chose to use PQ whether effective or not
    }
    return false;       // didn't use PQ
}


// ———  boring stuff: checkpointing and set-up  ———

void WHVivax::checkpoint(istream& stream){
    WHInterface::checkpoint(stream);
}
void WHVivax::checkpoint(ostream& stream){
    WHInterface::checkpoint(stream);
}

void WHVivax::init(){
    latentp = TimeStep(InputData().getModel().getParameters().getLatentp());
    initNHypnozoites();
}

}
}
