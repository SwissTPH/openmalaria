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

#include "WithinHost/CommonWithinHost.h"
#include "WithinHost/Diagnostic.h"
#include "WithinHost/Genotypes.h"
#include "util/errors.h"
#include "PopulationStats.h"
#include "util/AgeGroupInterpolation.h"
#include "util/random.h"
#include "util/StreamValidator.h"
#include "schema/scenario.h"

#include <boost/algorithm/string.hpp>

using namespace std;

namespace OM {
namespace WithinHost {

CommonInfection* (* CommonWithinHost::createInfection) (uint32_t protID);
CommonInfection* (* CommonWithinHost::checkpointedInfection) (istream& stream);

double hetMassMultStdDev = std::numeric_limits<double>::signaling_NaN();
double minHetMassMult = std::numeric_limits<double>::signaling_NaN();
util::AgeGroupInterpolator massByAge;

bool reportInfectedOrPatentInfected = false;
bool reportInfectionsByGenotype = false;

// Only required for a drug monitoring HACK and could be removed:
vector<string> drugMonCodes;


// -----  Initialization  -----

void CommonWithinHost::init( const scnXml::Scenario& scenario ){
    const scnXml::Human human = scenario.getModel().getHuman();
    if( !human.getWeight().present() ){
        // Technically this is needed by the PK/PD and Molineaux models
        throw util::xml_scenario_error( "model->human->weight element required by certain models" );
    }
    massByAge.set( human.getWeight().get(), "weight" );
    hetMassMultStdDev = human.getWeight().get().getMultStdDev();
    // hetWeightMult must be large enough that birth weight is at least 0.5 kg:
    minHetMassMult = 0.5 / massByAge.eval( 0.0 );
    
    const scnXml::Monitoring& mon = scenario.getMonitoring();
    if( mon.getDrugConcentration().present() ){
        boost::split( drugMonCodes,
                      mon.getDrugConcentration().get().getDrugCodes(),
                      boost::is_any_of("," ) );
    }
    
    reportInfectedOrPatentInfected = mon::isUsedM(mon::MHR_INFECTIONS) ||
        mon::isUsedM(mon::MHR_PATENT_INFECTIONS);
    reportInfectionsByGenotype = mon::isUsedM(mon::MHR_INFECTED_GENOTYPE) ||
        mon::isUsedM(mon::MHR_PATENT_GENOTYPE) ||
        mon::isUsedM(mon::MHF_LOG_DENSITY_GENOTYPE);
}

CommonWithinHost::CommonWithinHost( double comorbidityFactor ) :
        WHFalciparum( comorbidityFactor ), pkpdModel(PkPd::PkPdModel::createPkPdModel ())
{
    assert( sim::oneTS() == sim::fromDays(1) || sim::oneTS() == sim::fromDays(5) );
    
    // Sample a weight heterogeneity factor
#ifndef NDEBUG
    int counter = 0;
#endif
    do {
        hetMassMultiplier = util::random::gauss( 1.0, hetMassMultStdDev );
#ifndef NDEBUG
        assert( counter < 100 );        // too many resamples: resamples should rarely be needed...
        ++counter;
#endif
    } while( hetMassMultiplier < minHetMassMult );
}

CommonWithinHost::~CommonWithinHost() {
    delete pkpdModel;
    pkpdModel = 0;
    for( list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end(); ++inf ){
        delete *inf;
    }
    infections.clear();
}

// -----  Simple infection adders/removers  -----

void CommonWithinHost::clearInfections( Treatments::Stages stage ){
    for (std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end();) {
        if( stage == Treatments::BOTH ||
            (stage == Treatments::LIVER && !(*inf)->bloodStage()) ||
            (stage == Treatments::BLOOD && (*inf)->bloodStage())
        ){
            delete *inf;
            inf = infections.erase( inf );
        }else{
            ++inf;
        }
    }
    numInfs = infections.size();
}

// -----  interventions -----

void CommonWithinHost::treatPkPd(size_t schedule, size_t dosages, double age){
    double mass = massByAge.eval( age ) * hetMassMultiplier;
    pkpdModel->prescribe( schedule, dosages, age, mass );
}
void CommonWithinHost::clearImmunity() {
    for (std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end(); ++inf) {
        (*inf)->clearImmunity();
    }
    m_cumulative_h = 0.0;
    m_cumulative_Y_lag = 0.0;
}
void CommonWithinHost::importInfection(){
    PopulationStats::totalInfections += 1;
    if( numInfs < MAX_INFECTIONS ){
        PopulationStats::allowedInfections += 1;
        m_cumulative_h += 1;
        numInfs += 1;
        // This is a hook, used by interventions. The newly imported infections
        // should use initial frequencies to select genotypes.
        vector<double> weights( 0 );        // zero length: signal to use initial frequencies
        infections.push_back(createInfection(Genotypes::sampleGenotype(weights)));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}


// -----  Density calculations  -----

void CommonWithinHost::update(int nNewInfs, vector<double>& genotype_weights,
        double ageInYears, double bsvFactor, ofstream& drugMon)
{
    // Cache total density for infectiousness calculations
    int y_lag_i = sim::ts0().moduloSteps(y_lag_len);
    for( size_t g = 0; g < Genotypes::N(); ++g ) m_y_lag.at(y_lag_i, g) = 0.0;
    for( std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end(); ++inf ){
        m_y_lag.at( y_lag_i, (*inf)->genotype() ) += (*inf)->getDensity();
    }
    
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    PopulationStats::totalInfections += nNewInfs;
    nNewInfs=min(nNewInfs,MAX_INFECTIONS-numInfs);
    PopulationStats::allowedInfections += nNewInfs;
    numInfs += nNewInfs;
    assert( numInfs>=0 && numInfs<=MAX_INFECTIONS );
    for ( int i=0; i<nNewInfs; ++i ) {
        infections.push_back(createInfection (Genotypes::sampleGenotype(genotype_weights)));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
    
    updateImmuneStatus ();

    totalDensity = 0.0;
    timeStepMaxDensity = 0.0;
    
    // As in AJTMH p22, cumulative_h (X_h + 1) doesn't include infections added
    // this time-step and cumulative_Y only includes past densities.
    double cumulative_h=m_cumulative_h;
    double cumulative_Y=m_cumulative_Y;
    m_cumulative_h += nNewInfs;
    
    bool treatmentLiver = treatExpiryLiver > sim::ts0();
    bool treatmentBlood = treatExpiryBlood > sim::ts0();
    double survivalFactor_part = bsvFactor * _innateImmSurvFact;
    
    double body_mass = massByAge.eval( ageInYears ) * hetMassMultiplier;
    
    for( SimTime now = sim::ts0(), end = sim::ts0() + sim::oneTS(); now < end; now += sim::oneDay() ){
        // every day, medicate drugs, update each infection, then decay drugs
        pkpdModel->medicate( body_mass );
        
        double sumLogDens = 0.0;
        
        for (std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end();) {
            // Note: this is only one treatment model; there is also the PK/PD model
            bool expires = ((*inf)->bloodStage() ? treatmentBlood : treatmentLiver);
            
            if( !expires ){     /* no expiry due to simple treatment model; do update */
                double survivalFactor = survivalFactor_part *
                    (*inf)->immunitySurvivalFactor(ageInYears, cumulative_h, cumulative_Y) *
                    pkpdModel->getDrugFactor((*inf)->genotype());
                // update, may result in termination of infection:
                expires = (*inf)->update(survivalFactor, now, body_mass);
            }
            
            if( expires ){
                delete *inf;
                inf = infections.erase(inf);        // inf points to next infection now so don't increment with ++inf
                --numInfs;
            } else {
                double density = (*inf)->getDensity();
                totalDensity += density;
                timeStepMaxDensity = max(timeStepMaxDensity, density);
                m_cumulative_Y += density;
                if( density > 0 ){
                    // Base 10 logarithms are usually used; +1 because it avoids negatives in output while having very little affect on high densities
                    sumLogDens += log10(1.0 + density);
                }
                ++inf;
            }
        }
        pkpdModel->decayDrugs ();
        
        if( drugMon.is_open() && sim::intervNow() >= sim::zero() ){
            drugMon << now << '\t' << sumLogDens;
            map<string,double> concentrations;
            pkpdModel->getConcentrations( concentrations );
            foreach( string& drugCode, drugMonCodes ){
                drugMon << '\t' << concentrations[drugCode];
            }
            drugMon << endl;
        }
    }
    
    util::streamValidate(totalDensity);
    assert( (boost::math::isfinite)(totalDensity) );        // inf probably wouldn't be a problem but NaN would be
}

void CommonWithinHost::addProphylacticEffects(const vector<double>& pClearanceByTime) {
    // this should actually be easy; it just isn't needed yet
    throw util::unimplemented_exception( "prophylactic effects on 1-day time step" );
}


// -----  Summarize  -----

// Used in summarizeInfs.
vector<CommonInfection*> sortedInfs;
struct InfGenotypeSorter {
    bool operator() (CommonInfection* i, CommonInfection* j){
        return i->genotype() < j->genotype();
    }
} infGenotypeSorter;

void CommonWithinHost::summarizeInfs( const Host::Human& human )const{
    if( infections.size() == 0 ) return;        // nothing to report
    mon::reportMHI( mon::MHR_INFECTED_HOSTS, human, 1 );
    if( reportInfectedOrPatentInfected ){
        for (std::list<CommonInfection*>::const_iterator inf =
            infections.begin(); inf != infections.end(); ++inf) {
            uint32_t genotype = (*inf)->genotype();
            mon::reportMHGI( mon::MHR_INFECTIONS, human, genotype, 1 );
            if( Monitoring::Survey::diagnostic().isPositive( (*inf)->getDensity() ) ){
                mon::reportMHGI( mon::MHR_PATENT_INFECTIONS, human, genotype, 1 );
            }
        }
    }
    if( reportInfectionsByGenotype ){
        // Instead of storing nInfs and total density by genotype we sort
        // infections by genotype and report each in sequence.
        // We don't sort in place since that would affect random number sampling
        // order when updating, and the monitoring system should not in my
        // opinion affect outputs (since it would make testing harder).
        sortedInfs.assign( infections.begin(), infections.end() );
        sort( sortedInfs.begin(), sortedInfs.end(), infGenotypeSorter );
        vector<CommonInfection*>::const_iterator inf = sortedInfs.begin();
        while( inf != sortedInfs.end() ){
            uint32_t genotype = (*inf)->genotype();
            double dens = 0.0;
            do{     // at start: genotype is that of the current infection, dens is 0
                dens += (*inf)->getDensity();
                ++inf;
            }while( inf != sortedInfs.end() && (*inf)->genotype() == genotype );
            // we had at least one infection of this genotype
            mon::reportMHGI( mon::MHR_INFECTED_GENOTYPE, human, genotype, 1 );
            if( Monitoring::Survey::diagnostic().isPositive(dens) ){
                mon::reportMHGI( mon::MHR_PATENT_GENOTYPE, human, genotype, 1 );
                mon::reportMHGF( mon::MHF_LOG_DENSITY_GENOTYPE, human, genotype, log(dens) );
            }
        }
    }
}


void CommonWithinHost::checkpoint (istream& stream) {
    WHFalciparum::checkpoint (stream);
    hetMassMultiplier & stream;
    (*pkpdModel) & stream;
    for (int i = 0; i < numInfs; ++i) {
        infections.push_back (checkpointedInfection (stream));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}

void CommonWithinHost::checkpoint (ostream& stream) {
    WHFalciparum::checkpoint (stream);
    hetMassMultiplier & stream;
    (*pkpdModel) & stream;
    for (std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end(); ++inf) {
        (**inf) & stream;
    }
}
}
}
