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

#include "WithinHost/WHVivax.h"
#include "WithinHost/Pathogenesis/PathogenesisModel.h"
#include "WithinHost/Treatments.h"
#include "WithinHost/Genotypes.h"
#include "mon/reporting.h"
#include "util/random.h"
#include "util/errors.h"
#include "util/checkpoint_containers.h"
#include "util/timeConversions.h"
#include <schema/scenario.h>
#include <algorithm>
#include <limits>
#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;

// ———  parameters  ———

// Set from the parameters block:
SimTime latentP;       // attribute on parameters block

// Set from <vivax .../> element:
double probBloodStageInfectiousToMosq = numeric_limits<double>::signaling_NaN();
int maxNumberHypnozoites = -1;
double baseNumberHypnozoites = numeric_limits<double>::signaling_NaN();
int latentRelapseDays1stRelease = 0;
int latentRelapseDays2ndRelease = 0;
double pSecondRelease = numeric_limits<double>::signaling_NaN();
double muFirstHypnozoiteRelease = numeric_limits<double>::signaling_NaN();
double sigmaFirstHypnozoiteRelease = numeric_limits<double>::signaling_NaN();
double muSecondHypnozoiteRelease = numeric_limits<double>::signaling_NaN();
double sigmaSecondHypnozoiteRelease = numeric_limits<double>::signaling_NaN();
SimTime bloodStageProtectionLatency;
double bloodStageLengthWeibullScale = numeric_limits<double>::signaling_NaN();  // units: days
double bloodStageLengthWeibullShape = numeric_limits<double>::signaling_NaN();
double pPrimaryA = numeric_limits<double>::signaling_NaN(),
    pPrimaryB = numeric_limits<double>::signaling_NaN();
double pRelapseOneA = numeric_limits<double>::signaling_NaN(),
    pRelapseOneB = numeric_limits<double>::signaling_NaN();
double pRelapseTwoA = numeric_limits<double>::signaling_NaN(),
    pRelapseTwoB = numeric_limits<double>::signaling_NaN();
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
    assert(baseNumberHypnozoites <= 1 && baseNumberHypnozoites >= 0);
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
    if(isnan(pSecondRelease)) {
        pSecondRelease = 0.0;
    }
    bool isFirstRelease =
            (pSecondRelease == 0.0) ? true :
            ((pSecondRelease == 1.0) ? false :
            !random::bernoulli(pSecondRelease) );
    
    double mu, sigma;
    int latentRelapseDays;
    if (isFirstRelease){
        // only calculate a random delay from firstRelease distribution
        mu = muFirstHypnozoiteRelease;
        sigma = sigmaFirstHypnozoiteRelease;
        latentRelapseDays = latentRelapseDays1stRelease;
    } else {
        // only calculate a random delay from secondRelease distribution
        mu = muSecondHypnozoiteRelease;
        sigma = sigmaSecondHypnozoiteRelease;
        assert(!isnan(latentRelapseDays2ndRelease));
        latentRelapseDays = latentRelapseDays2ndRelease;
    }
    
    double liverStageMaximumDays = 16.0*30.0; // maximum of about 16 months in liver stage
    double delay = numeric_limits<double>::quiet_NaN();       // in days
    int count = 0;
    int maxcount = pow(10,6);
    
    do{
        delay = util::random::log_normal( mu, sigma );
        count += 1;
    }while( (delay > liverStageMaximumDays || delay < 0.0 ) && count < maxcount );
    if (count == maxcount) {
        throw util::xml_scenario_error( "<vivax><hypnozoiteRelease>  [random delay calculation causes probably an indefinite loop]:\n The hypnozoite release distribution seems off, sigma of secondRelease could be too high. We except the hypnozoite to reside a maximum of 16 months in the liver stage. Sigma choose well, dear padawan." );
    }
    assert( delay >= 0 && delay < liverStageMaximumDays );
    return sim::roundToTSFromDays( delay + latentRelapseDays );
}


// ———  per-brood code  ———

VivaxBrood::VivaxBrood( WHVivax *host ) :
        primaryHasStarted( false ),
        relapseHasStarted( false ),
        hadEvent( false ),
        hadRelapse( false )
{
    set<SimTime> releases;     // used to initialise releaseDates; a set is better to use now but a vector later
    
    // primary blood stage plus hypnozoites (relapses)
    releases.insert( sim::ts0() + latentP );
    int numberHypnozoites = sampleNHypnozoites();
    for( int i = 0; i < numberHypnozoites; ){
        SimTime randomReleaseDelay = sampleReleaseDelay();
        SimTime timeToRelease = sim::ts0() + latentP + randomReleaseDelay;
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
    relapseHasStarted & stream;
    hadEvent & stream;
    hadRelapse & stream;
}
VivaxBrood::VivaxBrood( istream& stream ){
    releaseDates & stream;
    bloodStageClearDate & stream;
    primaryHasStarted & stream;
    relapseHasStarted & stream;
    hadEvent & stream;
    hadRelapse & stream;
}


VivaxBrood::UpdResult VivaxBrood::update(){
    if( bloodStageClearDate == sim::ts0() ){
        //NOTE: this effectively means that both asexual and sexual stage
        // parasites self-terminate. It also means the immune system can
        // protect against new blood-stage infections for a short time.
    }
    
    UpdResult result;
    while( releaseDates.size() > 0 && releaseDates.back() == sim::ts0() ){
        releaseDates.pop_back();
        
#ifdef WHVivaxSamples
        if( sampleBrood == this ){
            cout << "Time\t" << sim::ts0();
            for( vector<SimTime>::const_iterator it = releaseDates.begin(); it != releaseDates.end(); ++it )
                cout << '\t' << *it;
            cout << endl;
        }
#endif
        
        // an existing or recently terminated blood stage from the same brood
        // protects against a newly released Hypnozoite
        //NOTE: this is an immunity effect: should there be no immunity when a blood stage first emerges?
        if( bloodStageClearDate + bloodStageProtectionLatency >= sim::ts0() ) continue;
        
        if( !relapseHasStarted && primaryHasStarted ){
            relapseHasStarted = true;
            result.newRelapseBS = true;
        }
        if( !primaryHasStarted ){
            primaryHasStarted = true;
            result.newPrimaryBS = true;
        }
        result.newBS = true;
        
        double lengthDays = random::weibull( bloodStageLengthWeibullScale, bloodStageLengthWeibullShape );
        bloodStageClearDate = sim::ts0() + sim::roundToTSFromDays( lengthDays );
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

WHVivax::WHVivax( double comorbidityFactor ) :
    cumPrimInf(0),
    pEvent( numeric_limits<double>::quiet_NaN() ),
    pFirstRelapseEvent( numeric_limits<double>::quiet_NaN() )
{
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

double WHVivax::probTransmissionToMosquito( double tbvFactor, double *sumX )const{
    assert( WithinHost::Genotypes::N() == 1 );
    assert( sumX == 0 );
    for (list<VivaxBrood>::const_iterator inf = infections.begin();
         inf != infections.end(); ++inf)
    {
        if( inf->isPatent() ){
            // we have gametocytes from at least one brood
            return probBloodStageInfectiousToMosq * tbvFactor;
        }
    }
    return 0;   // no gametocytes
}
double WHVivax::pTransGenotype(double pTrans, double sumX, size_t genotype){
    throw util::unimplemented_exception("genotype tracking for vivax");
}

bool WHVivax::summarize(const Host::Human& human) const{
    if( infections.size() == 0 ) return false;  // no infections: not patent, nothing to report
    mon::reportMHI( mon::MHR_INFECTED_HOSTS, human, 1 );
    bool patentHost = false;
    // (patent) infections are reported by genotype, even though we don't have
    // genotype in this model
    mon::reportMHGI( mon::MHR_INFECTIONS, human, 0, infections.size() );
    for (list<VivaxBrood>::const_iterator inf = infections.begin();
         inf != infections.end(); ++inf) 
    {
        if (inf->isPatent()){
            mon::reportMHGI( mon::MHR_PATENT_INFECTIONS, human, 0, 1 );
            patentHost = true;
        }
    }
    if( patentHost ) mon::reportMHI( mon::MHR_PATENT_HOSTS, human, 1 );
    return patentHost;
}

void WHVivax::importInfection(){
    // this means one new liver stage infection, which can result in multiple blood stages
    infections.push_back( VivaxBrood( this ) );
}

void WHVivax::update(int nNewInfs, vector<double>&,
        double ageInYears, double)
{
    // create new infections, letting the constructor do the initialisation work:
    for( int i = 0; i < nNewInfs; ++i )
        infections.push_back( VivaxBrood( this ) );
    
    // update infections
    // NOTE: currently no BSV model
    morbidity = Pathogenesis::NONE;
    uint32_t oldCumInf = cumPrimInf;
    double oldpEvent = ( isnan(pEvent))? 1.0 : pEvent;
    // always use the first relapse probability for following relapses as a factor
    double oldpRelapseEvent = ( isnan(pFirstRelapseEvent))? 1.0 : pFirstRelapseEvent;
    list<VivaxBrood>::iterator inf = infections.begin();
    while( inf != infections.end() ){
        VivaxBrood::UpdResult result = inf->update();
        if( result.newPrimaryBS ) cumPrimInf += 1;
        
        if( result.newBS ){
            // Sample for each new blood stage infection: the chance of some
            // clinical event.
            
            bool clinicalEvent;
            if( result.newPrimaryBS ){
                // Blood stage is primary. oldCumInf wasn't updated yet.
                double pPrimaryInfEvent = pPrimaryA * pPrimaryB / (pPrimaryB+oldCumInf);
                clinicalEvent = random::bernoulli( pPrimaryInfEvent );
                inf->setHadEvent( clinicalEvent );
            } else if ( result.newRelapseBS ){
                // Blood stage is a relapse. oldCumInf wasn't updated yet.
                double pFirstRelapseEvent = oldpEvent * (pRelapseOneA * pRelapseOneB / (pRelapseOneB + (oldCumInf-1)));
                clinicalEvent = random::bernoulli( pFirstRelapseEvent );
                inf->setHadRelapse( clinicalEvent );
            }else{
                // Subtract 1 from oldCumInf to not count the current brood in
                // the number of cumulative primary infections.
                if( inf->hasHadRelapse() ){
                    double pSecondRelapseEvent = oldpRelapseEvent * (pRelapseTwoA * pRelapseTwoB / (pRelapseTwoB + (oldCumInf-1)));
                    clinicalEvent = random::bernoulli( pSecondRelapseEvent );
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

bool WHVivax::diagnosticResult( const Diagnostic& diagnostic ) const{
    //TODO(monitoring): this shouldn't ignore the diagnostic (especially since
    // it should always return true if diagnostic.density=0)
    for (list<VivaxBrood>::const_iterator inf = infections.begin();
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

void WHVivax::treatment( Host::Human& human, TreatmentId treatId ){
    //TODO: something less ugly than this. Possibly we should revise code to
    // make Pf and Pv treatment models work more similarly.
    // For now we rely on the check in Treatments::Treatments(...).
    
    // This means clear blood stage infection(s) but not liver stage.
    for( list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
        it->treatmentBS();
    }
    
    // triggered intervention deployments:
    const Treatments& treat = Treatments::select( treatId );
    treat.deploy( human,
                  mon::Deploy::TREAT,
                  interventions::VaccineLimits(/*default initialise: no limits*/) );
}
bool WHVivax::optionalPqTreatment(){
    // PQ clears liver stages. We don't worry about the effect of PQ on
    // gametocytes, because these are always cleared by blood-stage drugs with
    // Vivax, and PQ is not given without BS drugs. NOTE: this ignores drug failure.
    if (pReceivePQ > 0.0 && !noPQ && random::bernoulli(pReceivePQ)){
        if( random::bernoulli(effectivenessPQ) ){
            for( list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
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
        for( list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
            it->treatmentLS();
        }
    }
    
    // there probably will be blood-stage treatment
    if( timeBlood < sim::zero() ){
        for( list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
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
        infections.push_back( VivaxBrood( stream ) );
    }
    noPQ & stream;
    int morbidity_i;
    morbidity_i & stream;
    morbidity = static_cast<Pathogenesis::State>( morbidity_i );
    cumPrimInf & stream;
    pEvent & stream;
    pFirstRelapseEvent & stream;
}
void WHVivax::checkpoint(ostream& stream){
    WHInterface::checkpoint(stream);
    infections.size() & stream;
    for( list<VivaxBrood>::iterator it = infections.begin(); it != infections.end(); ++it ){
        it->checkpoint( stream );
    }
    noPQ & stream;
    static_cast<int>( morbidity ) & stream;
    cumPrimInf & stream;
    pEvent & stream;
    pFirstRelapseEvent & stream;
}

char const*const not_impl = "feature not available in Vivax model";
void WHVivax::treatPkPd(size_t schedule, size_t dosages, double age){
    throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHVivax::getTotalDensity() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHVivax::getCumulative_h() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }
double WHVivax::getCumulative_Y() const{ throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }

void WHVivax::init( const OM::Parameters& parameters, const scnXml::Model& model ){
    try{
        //NOTE: if XSD is changed, this should not have a default unit
        latentP = UnitParse::readShortDuration( model.getParameters().getLatentp(), UnitParse::STEPS );
    }catch( const util::format_error& e ){
        throw util::xml_scenario_error( string("model/parameters/latentP: ").append(e.message()) );
    }
    if( !model.getVivax().present() )
        throw util::xml_scenario_error( "no vivax model description in scenario XML" );
    const scnXml::Vivax& elt = model.getVivax().get();
    probBloodStageInfectiousToMosq = elt.getProbBloodStageInfectiousToMosq().getValue();
    maxNumberHypnozoites = elt.getHypnozoiteRelease().getNumberHypnozoites().getMax();
    baseNumberHypnozoites = elt.getHypnozoiteRelease().getNumberHypnozoites().getBase();
    const scnXml::HypnozoiteRelease& hr = elt.getHypnozoiteRelease();
    latentRelapseDays1stRelease = hr.getFirstRelease().getLatentRelapseDays();
    muFirstHypnozoiteRelease = hr.getFirstRelease().getMu();
    sigmaFirstHypnozoiteRelease = hr.getFirstRelease().getSigma();
    if(hr.getSecondRelease().present()){
        latentRelapseDays2ndRelease = hr.getSecondRelease().get().getLatentRelapseDays();
        muSecondHypnozoiteRelease = hr.getSecondRelease().get().getMu();
        sigmaSecondHypnozoiteRelease = hr.getSecondRelease().get().getSigma();
        pSecondRelease = hr.getPSecondRelease();
        if( pSecondRelease == 0.0){
            std::cerr << "Warning: probability of second release is set to zero, although secondRelease element is present, will only calculate for a first release." << endl;
        }
        assert( pSecondRelease >= 0 && pSecondRelease <= 1 );
    }
    bloodStageProtectionLatency = sim::roundToTSFromDays( elt.getBloodStageProtectionLatency().getValue() );
    bloodStageLengthWeibullScale = elt.getBloodStageLengthDays().getWeibullScale();
    bloodStageLengthWeibullShape = elt.getBloodStageLengthDays().getWeibullShape();
    
    const scnXml::ClinicalEvents& ce = elt.getClinicalEvents();
    pPrimaryA = ce.getPPrimaryInfection().getA();
    pPrimaryB = ce.getPPrimaryInfection().getB();
    pRelapseOneA = ce.getPRelapseOne().getA();
    pRelapseOneB = ce.getPRelapseOne().getB();
    pRelapseTwoA = ce.getPRelapseTwoPlus().getA();
    pRelapseTwoB = ce.getPRelapseTwoPlus().getB();
    pEventIsSevere = ce.getPEventIsSevere().getValue();
    
    initNHypnozoites();
    Pathogenesis::PathogenesisModel::init( parameters, model.getClinical(), true );
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
