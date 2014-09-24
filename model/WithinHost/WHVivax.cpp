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
#include "WithinHost/Pathogenesis/PathogenesisModel.h"
#include "WithinHost/Treatments.h"
#include "util/random.h"
#include "util/errors.h"
#include <schema/scenario.h>
#include <algorithm>
#include <limits>
#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;
using namespace Monitoring;
using boost::ptr_list;

// ———  parameters  ———

// Set from the parameters block:
SimTime latentP;       // attribute on parameters block

// Set from <vivax .../> element:
double probBloodStageInfectiousToMosq = numeric_limits<double>::signaling_NaN();
int maxNumberHypnozoites = -1;
double baseNumberHypnozoites = numeric_limits<double>::signaling_NaN();
double muReleaseHypnozoite = numeric_limits<double>::signaling_NaN();   // units: days
double sigmaReleaseHypnozoite = numeric_limits<double>::signaling_NaN();
double minReleaseHypnozoite;    // units: days
SimTime bloodStageProtectionLatency;
double bloodStageLengthWeibullScale = numeric_limits<double>::signaling_NaN();  // units: days
double bloodStageLengthWeibullShape = numeric_limits<double>::signaling_NaN();
double pEventPrimA = numeric_limits<double>::signaling_NaN(),
    pEventPrimB = numeric_limits<double>::signaling_NaN();
double pEventSecA = numeric_limits<double>::signaling_NaN(),
    pEventSecB = numeric_limits<double>::signaling_NaN();
double pEventIsSevere = numeric_limits<double>::signaling_NaN();

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
SimTime sampleReleaseDelay(){
    double delay;       // in days
    do{
        delay = util::random::log_normal( muReleaseHypnozoite, sigmaReleaseHypnozoite );
    }while( delay < minReleaseHypnozoite );
    return sim::roundToTSFromDays( delay );
}


// ———  per-brood code  ———

VivaxBrood::VivaxBrood( WHVivax *host ) :
        primaryHasStarted( false ),
        hadEvent( false )
{
    set<SimTime> releases;     // used to initialise releaseDates; a set is better to use now but a vector later
    
    // primary blood stage plus hypnozoites (relapses)
    releases.insert( sim::now1() + latentP );
    int numberHypnozoites = sampleNHypnozoites();
    for( int i = 0; i < numberHypnozoites; ){
        SimTime timeToRelease = sim::now1() + latentP + sampleReleaseDelay();
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
        for( vector<SimTime>::const_iterator it = releaseDates.begin(); it != releaseDates.end(); ++it )
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
    hadEvent & stream;
}
VivaxBrood::VivaxBrood( istream& stream ){
    releaseDates & stream;
    bloodStageClearDate & stream;
    primaryHasStarted & stream;
    hadEvent & stream;
}


VivaxBrood::UpdResult VivaxBrood::update(){
    if( bloodStageClearDate == sim::now1() ){
        //NOTE: this effectively means that both asexual and sexual stage
        // parasites self-terminate. It also means the immune system can
        // protect against new blood-stage infections for a short time.
    }
    
    UpdResult result;
    while( releaseDates.size() > 0 && releaseDates.back() == sim::now1() ){
        releaseDates.pop_back();
        
#ifdef WHVivaxSamples
        if( sampleBrood == this ){
            cout << "Time\t" << sim::now1();
            for( vector<SimTime>::const_iterator it = releaseDates.begin(); it != releaseDates.end(); ++it )
                cout << '\t' << *it;
            cout << endl;
        }
#endif
        
        // an existing or recently terminated blood stage from the same brood
        // protects against a newly released Hypnozoite
        //NOTE: this is an immunity effect: should there be no immunity when a blood stage first emerges?
        if( bloodStageClearDate + bloodStageProtectionLatency
            >= sim::now1() ) continue;
        
        if( !primaryHasStarted ){
            primaryHasStarted = true;
            result.newPrimaryBS = true;
        }
        result.newBS = true;
        
        double lengthDays = random::weibull( bloodStageLengthWeibullScale, bloodStageLengthWeibullShape );
        bloodStageClearDate = sim::now1() + sim::roundToTSFromDays( lengthDays );
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
    bloodStageClearDate = sim::never();
}

void VivaxBrood::treatmentLS(){
    releaseDates.clear();       // 100% clearance
    
    /* partial clearance code, in case of need:
    vector<SimTime> survivingZoites;
    survivingZoites.reserve( releaseDates.size() );   // maximum size we need
    for( vector<SimTime>::const_iterator it = releaseDates.begin(); it != releaseDates.end(); ++it ){
        if( !random::bernoulli( pClearEachHypnozoite ) ){
            survivingZoites.push_back( *it );    // copy            
        }
    }
    swap( releaseDates, survivingZoites );
    */
}


// ———  per-host code  ———

WHVivax::WHVivax( double comorbidityFactor ) : cumPrimInf(0) {
    if( comorbidityFactor != 1.0 )
#ifdef WHVivaxSamples
    if( sampleHost == 0 ){
        sampleHost = this;
        cout << "New host" << endl;
    }
#endif
        throw TRACED_EXCEPTION_DEFAULT( "This vivax model cannot be used with comorbidity heterogeneity" );
    noPQ = ( pHetNoPQ > 0.0 && random::bernoulli(pHetNoPQ) );
}

WHVivax::~WHVivax(){
#ifdef WHVivaxSamples
    if( this == sampleHost ){
        sampleHost = 0;
        cout << "Host terminates" << endl;
    }
#endif
}

double WHVivax::probTransmissionToMosquito( double tbvFactor )const{
    for (ptr_list<VivaxBrood>::const_iterator inf = infections.begin();
         inf != infections.end(); ++inf)
    {
        if( inf->isPatent() ){
            // we have gametocytes from at least one brood
            return probBloodStageInfectiousToMosq * tbvFactor;
        }
    }
    return 0;   // no gametocytes
}

bool WHVivax::summarize(const Host::Human& human) {
    Survey& survey = Survey::current();
    InfectionCount count = countInfections();
    if( count.total != 0 ){
        survey
            .addInt( Report::MI_INFECTED_HOSTS, human, 1 )
            .addInt( Report::MI_INFECTIONS, human, count.total )
            .addInt( Report::MI_PATENT_INFECTIONS, human, count.patent );
    }
    if( count.patent > 0 ){
        survey
            .addInt( Report::MI_PATENT_HOSTS, human, 1 );
        return true;
    }
    return false;
}

void WHVivax::importInfection(){
    // this means one new liver stage infection, which can result in multiple blood stages
    infections.push_back( new VivaxBrood( this ) );
}

void WHVivax::update(int nNewInfs, double ageInYears, double, ofstream& drugMon){
    // create new infections, letting the constructor do the initialisation work:
    for( int i = 0; i < nNewInfs; ++i )
        infections.push_back( new VivaxBrood( this ) );
    
    // update infections
    // NOTE: currently no BSV model
    morbidity = Pathogenesis::NONE;
    uint32_t oldCumInf = cumPrimInf;
    ptr_list<VivaxBrood>::iterator inf = infections.begin();
    while( inf != infections.end() ){
        VivaxBrood::UpdResult result = inf->update();
        if( result.newPrimaryBS ) cumPrimInf += 1;
        
        if( result.newBS ){
            // Sample for each new blood stage infection: the chance of some
            // clinical event.
            
            bool clinicalEvent;
            if( result.newPrimaryBS ){
                // Blood stage is primary. oldCumInf wasn't updated yet.
                double pEvent = pEventPrimA * pEventPrimB / (pEventPrimB+oldCumInf);
                clinicalEvent = random::bernoulli( pEvent );
                inf->setHadEvent( clinicalEvent );
            }else{
                // Subtract 1 from oldCumInf to not count the current brood in
                // the number of cumulative primary infections.
                if( inf->hasHadEvent() ){
                    double pEvent = pEventSecA * pEventSecB / (pEventSecB + (oldCumInf-1));
                    clinicalEvent = random::bernoulli( pEvent );
                }else{
                    // If the primary infection did not cause an event, there
                    // is 0 chance of a secondary causing an event in our model.
                    clinicalEvent = false;
                }
            }
            
            if( clinicalEvent ){
                if( random::bernoulli( pEventIsSevere ) )
                    morbidity = static_cast<Pathogenesis::State>( morbidity | Pathogenesis::STATE_SEVERE );
                else
                    morbidity = static_cast<Pathogenesis::State>( morbidity | Pathogenesis::STATE_MALARIA );
            }
        }
        
        if( result.isFinished ) inf = infections.erase( inf );
        else ++inf;
    }
    
    //NOTE: currently we don't model co-infection or indirect deaths
    if( morbidity == Pathogenesis::NONE ){
        morbidity = Pathogenesis::PathogenesisModel::sampleNMF( ageInYears );
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

Pathogenesis::StatePair WHVivax::determineMorbidity(double ageYears){
    Pathogenesis::StatePair result;     // no indirect mortality in the vivax model
    result.state = morbidity;
    return result;
}

void WHVivax::clearImmunity(){
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
void WHVivax::treatment( Host::Human& human, TreatmentId treatId ){
    //TODO: something less ugly than this. Possibly we should revise code to
    // make Pf and Pv treatment models work more similarly.
    // For now we rely on the check in Treatments::Treatments(...).
    
    // This means clear blood stage infection(s) but not liver stage.
    for( ptr_list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
        it->treatmentBS();
    }
    
    // triggered intervention deployments:
    const Treatments& treat = Treatments::select( treatId );
    treat.deploy( human,
                  interventions::Deployment::TREAT,
                  interventions::VaccineLimits(/*default initialise: no limits*/) );
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
void WHVivax::treatSimple(SimTime timeLiver, SimTime timeBlood){
    //TODO: this should be implemented properly (allowing effects on next
    // update instead of now)
    
    // liver-stage treatment is only via "Primaquine" option, if at all
    if( timeLiver != sim::zero() ){
        if( pReceivePQ > 0.0 ){
            throw util::xml_scenario_error("simple treatment for vivax liver "
            "stages is incompatible with case-management Primaquine option");
        }
        if( timeLiver >= sim::zero() )
            throw util::unimplemented_exception("simple treatment for vivax, except with timesteps=-1");
        for( ptr_list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
            it->treatmentLS();
        }
    }
    
    // there probably will be blood-stage treatment
    if( timeBlood < sim::zero() ){
        for( ptr_list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
            it->treatmentBS();
        }
    }else{
        if( timeBlood != sim::zero() )
            throw util::unimplemented_exception("simple treatment for vivax, except with timesteps=-1");
    }
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

char const*const not_impl = "feature not available in Vivax model";
void WHVivax::treatPkPd(size_t schedule, size_t dosages, double age){
    throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHVivax::getTotalDensity() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHVivax::getCumulative_h() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHVivax::getCumulative_Y() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }

void WHVivax::init( const OM::Parameters& parameters, const scnXml::Scenario& scenario ){
    //FIXME(schema): should be entered in days
    latentP = sim::fromTS(  scenario.getModel().getParameters().getLatentp() );
    if( !scenario.getModel().getVivax().present() )
        throw util::xml_scenario_error( "no vivax model description in scenario XML" );
    const scnXml::Vivax& elt = scenario.getModel().getVivax().get();
    probBloodStageInfectiousToMosq = elt.getProbBloodStageInfectiousToMosq().getValue();
    maxNumberHypnozoites = elt.getNumberHypnozoites().getMax();
    baseNumberHypnozoites = elt.getNumberHypnozoites().getBase();
    muReleaseHypnozoite = elt.getHypnozoiteReleaseDelayDays().getMu();
    sigmaReleaseHypnozoite = elt.getHypnozoiteReleaseDelayDays().getSigma();
    minReleaseHypnozoite = elt.getHypnozoiteReleaseDelayDays().getMin();
    bloodStageProtectionLatency = sim::roundToTSFromDays( elt.getBloodStageProtectionLatency().getValue() );
    bloodStageLengthWeibullScale = elt.getBloodStageLengthDays().getWeibullScale();
    bloodStageLengthWeibullShape = elt.getBloodStageLengthDays().getWeibullShape();
    
    pEventPrimA = elt.getPEventPrimary().getA();
    pEventPrimB = elt.getPEventPrimary().getB();
    pEventSecA = elt.getPEventSecondary().getA();
    pEventSecB = elt.getPEventSecondary().getB();
    pEventIsSevere = elt.getPEventIsSevere().getValue();
    
    initNHypnozoites();
    Pathogenesis::PathogenesisModel::init( parameters, scenario.getModel().getClinical(), true );
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
