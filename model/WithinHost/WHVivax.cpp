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
#include <algorithm>
#include <limits>

namespace OM {
namespace WithinHost {

using namespace OM::util;

//TODO: these parameters need to be initialised from XML, and possibly made
// static members of VivaxInfection / WHVivax
// Note: current values need replacing
double meanNumberHypnozoites = 5;       // have no idea what value to use
double meanTimeToHypReleaseInDays = 4;
double meanDurationBloodStageDays = 8;
double shapeParameterDurationBloodStage = 2;
double probPatentInfectMosquitoes = 1;
double probLiverStageTreatment = 0;
double efficacyLiverStageTreatment = 1;


struct TimeStepCompReverse{
    bool operator()( const TimeStep& a, const TimeStep& b ){
        return a < b;
    }
} timeStepCompR;
VivaxBrood::VivaxBrood(){
    int numberHypnozoites = random::poisson( meanNumberHypnozoites );
    for( int i = 0; i < numberHypnozoites; ++i ){
        double timeToReleaseDays = random::exponential( meanTimeToHypReleaseInDays );
        TimeStep timeToRelease = TimeStep::fromDaysNearest( timeToReleaseDays );
        hypReleaseTimes.push_back( TimeStep::simulation + timeToRelease );
    }
    
    // Sort by time to release, smallest (soonest) first. Explanation of code:
    // http://www.cplusplus.com/reference/algorithm/sort/
    sort( hypReleaseTimes.begin(), hypReleaseTimes.end(), timeStepCompR );
}

bool VivaxBrood::update( bool& anyNewBloodStage ){
    if( bloodStageClearTime == TimeStep::simulation )
        bloodStageClearTime = TimeStep::never;
    
    while( hypReleaseTimes.back() == TimeStep::simulation ){
        hypReleaseTimes.pop_back();
        // note that hypnozoite release during an existing blood stage infection does nothing
        if( isPatent() ) continue;
        
        anyNewBloodStage = true;
        double lengthDays = random::weibull( meanDurationBloodStageDays, shapeParameterDurationBloodStage );
        bloodStageClearTime = TimeStep::simulation + TimeStep::fromDaysNearest( lengthDays );
    }
    
    bool isFinished = (!isPatent()) && (hypReleaseTimes.size() == 0);
    return isFinished;
}

void VivaxBrood::treatment( double pClearEachHypnozoite ){
    bloodStageClearTime = TimeStep::never;      // clear blood stage
    
    if( pClearEachHypnozoite <= 0 )
        return;     // no liver stage effect
    if( pClearEachHypnozoite >= 1 ){
        hypReleaseTimes.clear();        // full liver stage clearance
        return;
    }
    
    // partial liver stage clearance:
    vector<TimeStep> newVec;
    newVec.reserve( hypReleaseTimes.size() );   // maximum size we need
    for( vector<TimeStep>::const_iterator it = hypReleaseTimes.begin(); it != hypReleaseTimes.end(); ++it ){
        if( !random::bernoulli( pClearEachHypnozoite ) ){
            newVec.push_back( *it );    // copy            
        }
    }
    swap( hypReleaseTimes, newVec );
}

void WHVivax::init(){
    //TODO: read params from XML
}

WHVivax::WHVivax( double comorbidityFactor ){
    if( comorbidityFactor != 1.0 )
        throw TRACED_EXCEPTION_DEFAULT( "This vivax model cannot be used with comorbidity heterogeneity" );
    //TODO: do we need to do anything here?
    throw xml_scenario_error( "vivax model code is not yet functional" );
}

WHVivax::~WHVivax(){

}

double WHVivax::probTransmissionToMosquito(TimeStep ageTimeSteps, double tbvFactor) const{
    if( !diagnosticDefault() ) return 0;        // no blood stage infections: no transmission
    
    // Otherwise we have at least one blood stage. We use a simple model and
    // ignore lag, the effect of Primaquine, exact densities, and recombination.
    return probPatentInfectMosquitoes * tbvFactor;
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

void WHVivax::update(int nNewInfs, double ageInYears, double){
    // create new infections, letting the constructor do the initialisation work:
    infections.resize( infections.size() + nNewInfs );
    
    // update infections
    // TODO: how do blood stage vaccines work?
    bool anyNewBloodStage = false;
    std::list<VivaxBrood>::iterator inf = infections.begin();
    while( inf != infections.end() ){
        bool isFinished = inf->update( anyNewBloodStage );
        if( isFinished ) inf = infections.erase( inf );
        else ++inf;
    }
    
    morbidity = Pathogenesis::NONE;
    //TODO: does each hypnozoite release have a chance of causing morbidity or
    // does this happen just once per timestep when at least one hypnozoite releases?
    if( anyNewBloodStage ){
        //TODO: model?
        double xyz = numeric_limits<double>::quiet_NaN();
        if( random::bernoulli( xyz ) )
            morbidity = Pathogenesis::SICK;  //TODO: or severe?
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

Pathogenesis::State WHVivax::determineMorbidity(double ageYears){
    return morbidity;
}

void WHVivax::clearImmunity(){
    //TODO: reset all immunity (what do we have?)
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
    // This means clear blood stage infection(s), and
    // possibly use Primaquine to clear some hypnozoites from the liver.
    
    double liverStageEfficacy = 0;
    if( random::bernoulli( probLiverStageTreatment ) ){
        liverStageEfficacy = efficacyLiverStageTreatment;
    }
    
    for( list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
        it->treatment( liverStageEfficacy );
    }
}

void WHVivax::checkpoint(istream& stream){
    WHInterface::checkpoint(stream);
}
void WHVivax::checkpoint(ostream& stream){
    WHInterface::checkpoint(stream);
}

}
}
