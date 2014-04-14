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
#include <schema/healthSystem.h>
#include <algorithm>
#include <limits>
#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;
using boost::ptr_list;

// ———  parameters  ———

// Set from the parameters block:
TimeStep latentp;       // attribute on parameters block

// Set from <vivax .../> element:
double probBloodStageInfectiousToMosq = numeric_limits<double>::signaling_NaN();
int maxNumberHypnozoites = -1;
double baseNumberHypnozoites = numeric_limits<double>::signaling_NaN();
double muReleaseHypnozoite = numeric_limits<double>::signaling_NaN();   // units: days
double sigmaReleaseHypnozoite = numeric_limits<double>::signaling_NaN();
double minReleaseHypnozoite;    // units: days
TimeStep bloodStageProtectionLatency;
double bloodStageLengthWeibullScale = numeric_limits<double>::signaling_NaN();  // units: days
double bloodStageLengthWeibullShape = numeric_limits<double>::signaling_NaN();

// Set from healthSystem element:
double pHetNoPQ = numeric_limits<double>::signaling_NaN();
double pReceivePQ = numeric_limits<double>::signaling_NaN();
double effectivenessPQ = numeric_limits<double>::signaling_NaN();

#ifdef WHVivaxSamples
WHVivax *sampleHost = 0;
VivaxBrood *sampleBrood = 0;
#endif


// ———  individual models  ———

// number of hypnozoites per brood:
map<double,int> nHypnozoitesProbMap;
void initNHypnozoites(){
    double total = 0.0;
    for( int n = 0; n <= maxNumberHypnozoites; ++n )
        total += pow( baseNumberHypnozoites, n );
    
    double cumP = 0.0;
    for( int n = 0; n <= maxNumberHypnozoites; ++n ){
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
    double delay;       // in days
    do{
        delay = util::random::log_normal( muReleaseHypnozoite, sigmaReleaseHypnozoite );
    }while( delay < minReleaseHypnozoite );
    return TimeStep::fromDaysNearest( delay );
}


// ———  per-brood code  ———

VivaxBrood::VivaxBrood( WHVivax *host ) :
        primaryHasStarted( false )
{
    set<TimeStep> releases;     // used to initialise releaseDates; a set is better to use now but a vector later
    
    // primary blood stage plus hypnozoites (relapses)
    releases.insert( TimeStep::simulation + latentp );
    int numberHypnozoites = sampleNHypnozoites();
    for( int i = 0; i < numberHypnozoites; ){
        TimeStep timeToRelease = TimeStep::simulation + latentp + sampleReleaseDelay();
        bool inserted = releases.insert( timeToRelease ).second;
        if( inserted ) ++i;     // successful
        // else: sample clash with an existing release date, so resample
    }
    
    // Copy times to the vector, backwards (smallest last):
    releaseDates.insert( releaseDates.end(), releases.rbegin(), releases.rend() );
    
#ifdef WHVivaxSamples
    if( sampleHost == host && sampleBrood == 0 ){
        sampleBrood = this;
        cout << "New sample brood";
        for( vector<TimeStep>::const_iterator it = releaseDates.begin(); it != releaseDates.end(); ++it )
            cout << '\t' << *it;
        cout << endl;
    }
#endif
}
VivaxBrood::~VivaxBrood(){
#ifdef WHVivaxSamples
    if( sampleBrood == this ){
        sampleBrood = 0;
        cout << "Brood terminated" << endl;
    }
#endif
}

void VivaxBrood::checkpoint( ostream& stream ){
    releaseDates & stream;
    bloodStageClearDate & stream;
    primaryHasStarted & stream;
}
VivaxBrood::VivaxBrood( istream& stream ){
    releaseDates & stream;
    bloodStageClearDate & stream;
    primaryHasStarted & stream;
}


VivaxBrood::UpdatePair VivaxBrood::update( bool& anyNewBloodStage ){
    if( bloodStageClearDate == TimeStep::simulation ){
        //NOTE: this effectively means that both asexual and sexual stage
        // parasites self-terminate. It also means the immune system can
        // protect against new blood-stage infections for a short time.
    }
    
    UpdatePair result;
    result.newPrimaryBS = false;
    while( releaseDates.size() > 0 && releaseDates.back() == TimeStep::simulation ){
        releaseDates.pop_back();
        if( !primaryHasStarted ){
            primaryHasStarted = true;
            result.newPrimaryBS = true;
        }
        
#ifdef WHVivaxSamples
        if( sampleBrood == this ){
            cout << "Time\t" << TimeStep::simulation;
            for( vector<TimeStep>::const_iterator it = releaseDates.begin(); it != releaseDates.end(); ++it )
                cout << '\t' << *it;
            cout << endl;
        }
#endif
        
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
    
    result.isFinished = (!isPatent()) && (releaseDates.size() == 0);
    return result;
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

WHVivax::WHVivax() : cumPrimInf(0) {
    noPQ = ( pHetNoPQ > 0.0 && random::bernoulli(pHetNoPQ) );
#ifdef WHVivaxSamples
    if( sampleHost == 0 ){
        sampleHost = this;
        cout << "New host" << endl;
    }
#endif
}

void WHVivax::setComorbidityFactor(double factor){
    if( factor != 1.0 )
        throw TRACED_EXCEPTION_DEFAULT( "This vivax model cannot be used with comorbidity heterogeneity" );
}

WHVivax::~WHVivax(){
#ifdef WHVivaxSamples
    if( this == sampleHost ){
        sampleHost = 0;
        cout << "Host terminates" << endl;
    }
#endif
}

double WHVivax::probTransmissionToMosquito(TimeStep ageTimeSteps, double tbvEfficacy) const{
    for (ptr_list<VivaxBrood>::const_iterator inf = infections.begin();
         inf != infections.end(); ++inf)
    {
        if( inf->isPatent() ){
            // we have gametocytes from at least one brood
            return probBloodStageInfectiousToMosq * (1.0 - tbvEfficacy );
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
    infections.push_back( new VivaxBrood( this ) );
}

void WHVivax::update(int nNewInfs, double ageInYears, double BSVEfficacy){
    // create new infections, letting the constructor do the initialisation work:
    for( int i = 0; i < nNewInfs; ++i )
        infections.push_back( new VivaxBrood( this ) );
    
    // update infections
    // NOTE: currently no BSV model
    bool anyNewBloodStage = false;
    ptr_list<VivaxBrood>::iterator inf = infections.begin();
    while( inf != infections.end() ){
        VivaxBrood::UpdatePair result = inf->update( anyNewBloodStage );
        if( result.newPrimaryBS ) cumPrimInf += 1;
        if( result.isFinished ) inf = infections.erase( inf );
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
    for (ptr_list<VivaxBrood>::const_iterator inf = infections.begin();
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
    for (ptr_list<VivaxBrood>::const_iterator inf = infections.begin(); inf != infections.end(); ++inf) {
        if (inf->isPatent())
            count.patent += 1;
    }
    return count;
}
void WHVivax::effectiveTreatment(){
    // This means clear blood stage infection(s) but not liver stage.
    for( ptr_list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
        it->treatmentBS();
    }
}
bool WHVivax::optionalPqTreatment(){
    // PQ clears liver stages. We don't worry about the effect of PQ on
    // gametocytes, because these are always cleared by blood-stage drugs with
    // Vivax, and PQ is not given without BS drugs. NOTE: this ignores drug failure.
    if (pReceivePQ > 0.0 && !noPQ && random::bernoulli(pReceivePQ)){
        if( random::bernoulli(effectivenessPQ) ){
            for( ptr_list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
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
    size_t len;
    len & stream;
    for( size_t i = 0; i < len; ++i ){
        infections.push_back( new VivaxBrood( stream ) );
    }
    noPQ & stream;
    int morbidity_i;
    morbidity_i & stream;
    morbidity = static_cast<Pathogenesis::State>( morbidity_i );
    cumPrimInf & stream;
}
void WHVivax::checkpoint(ostream& stream){
    WHInterface::checkpoint(stream);
    infections.size() & stream;
    for( ptr_list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
        it->checkpoint( stream );
    }
    noPQ & stream;
    static_cast<int>( morbidity ) & stream;
    cumPrimInf & stream;
}

void WHVivax::init(){
    latentp = TimeStep(InputData().getModel().getParameters().getLatentp());
    if( !InputData().getModel().getVivax().present() )
        throw util::xml_scenario_error( "no vivax model description in scenario XML" );
    const scnXml::Vivax& elt = InputData().getModel().getVivax().get();
    probBloodStageInfectiousToMosq = elt.getProbBloodStageInfectiousToMosq().getValue();
    maxNumberHypnozoites = elt.getNumberHypnozoites().getMax();
    baseNumberHypnozoites = elt.getNumberHypnozoites().getBase();
    muReleaseHypnozoite = elt.getHypnozoiteReleaseDelayDays().getMu();
    sigmaReleaseHypnozoite = elt.getHypnozoiteReleaseDelayDays().getSigma();
    minReleaseHypnozoite = elt.getHypnozoiteReleaseDelayDays().getMin();
    bloodStageProtectionLatency = TimeStep::fromDaysNearest( elt.getBloodStageProtectionLatency().getValue() );
    bloodStageLengthWeibullScale = elt.getBloodStageLengthDays().getWeibullScale();
    bloodStageLengthWeibullShape = elt.getBloodStageLengthDays().getWeibullShape();
    
    initNHypnozoites();
}
void WHVivax::setHSParameters(const scnXml::Primaquine& elt){
    if( pHetNoPQ != pHetNoPQ )// not set yet
        pHetNoPQ = elt.getPHumanCannotReceive().getValue();
    else{
        if( pHetNoPQ != elt.getPHumanCannotReceive().getValue() )
            throw util::xml_scenario_error( "changeHS cannot change pHumanCannotReceive value" );
    }
    pReceivePQ = elt.getPUseUncomplicated().getValue();
    effectivenessPQ = elt.getEffectivenessOnUse().getValue();
}

}
}
