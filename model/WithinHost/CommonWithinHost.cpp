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
#include "util/errors.h"
#include "PopulationStats.h"
#include "util/StreamValidator.h"

#include <boost/algorithm/string.hpp>

using namespace std;

namespace OM {
namespace WithinHost {

CommonInfection* (* CommonWithinHost::createInfection) (uint32_t protID);
CommonInfection* (* CommonWithinHost::checkpointedInfection) (istream& stream);

vector<string> drugMonCodes;


// -----  Initialization  -----

void CommonWithinHost::init(const scnXml::DrugConcentration& elt){
    boost::split( drugMonCodes, elt.getDrugCodes(), boost::is_any_of("," ) );
}

CommonWithinHost::CommonWithinHost( double comorbidityFactor ) :
        WHFalciparum( comorbidityFactor ), pkpdModel(PkPd::PkPdModel::createPkPdModel ())
{
    assert( sim::oneTS() == sim::fromDays(1) || sim::oneTS() == sim::fromDays(5) );
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
    pkpdModel->prescribe( schedule, dosages, age );
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
        infections.push_back(createInfection (pkpdModel->new_proteome_ID ()));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}


// -----  Density calculations  -----

void CommonWithinHost::update(int nNewInfs, double ageInYears, double bsvFactor, ofstream& drugMon) {
    // Cache total density for infectiousness calculations
    m_y_lag[sim::nowStepsMod(y_lag_len)] = totalDensity;
    
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    PopulationStats::totalInfections += nNewInfs;
    nNewInfs=min(nNewInfs,MAX_INFECTIONS-numInfs);
    PopulationStats::allowedInfections += nNewInfs;
    numInfs += nNewInfs;
    assert( numInfs>=0 && numInfs<=MAX_INFECTIONS );
    for ( int i=0; i<nNewInfs; ++i ) {
        infections.push_back(createInfection (pkpdModel->new_proteome_ID ()));
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

    bool treatmentLiver = treatExpiryLiver >= sim::now();
    bool treatmentBlood = treatExpiryBlood >= sim::now();
    double survivalFactor_part = bsvFactor * _innateImmSurvFact;
    
    for( SimTime now = sim::now(), end = sim::now() + sim::oneTS(); now < end; now += sim::oneDay() ){
        // every day, medicate drugs, update each infection, then decay drugs
        pkpdModel->medicate( ageInYears );
        
        double sumLogDens = 0.0;
        
        for (std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end();) {
            // Note: this is only one treatment model; there is also the PK/PD model
            bool expires = ((*inf)->bloodStage() ? treatmentBlood : treatmentLiver);
            
            if( !expires ){     /* no expiry due to simple treatment model; do update */
                double survivalFactor = survivalFactor_part *
                    (*inf)->immunitySurvivalFactor(ageInYears, cumulative_h, cumulative_Y) *
                    pkpdModel->getDrugFactor((*inf)->get_proteome_ID());
                // update, may result in termination of infection:
                expires = (*inf)->update(survivalFactor, now);
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
                    //NOTE: this is a provisional output. It should not, however, add log(0)=-inf to the sum.
                    sumLogDens += log(density);
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

WHInterface::InfectionCount CommonWithinHost::countInfections () const{
    InfectionCount count;       // constructor initialises counts to 0
    count.total = infections.size();
    for (std::list<CommonInfection*>::const_iterator inf = infections.begin(); inf != infections.end(); ++inf) {
        if (Diagnostic::default_.isPositive( (*inf)->getDensity() ) )
            count.patent += 1;
    }
    return count;
}


void CommonWithinHost::checkpoint (istream& stream) {
    WHFalciparum::checkpoint (stream);
    (*pkpdModel) & stream;
    for (int i = 0; i < numInfs; ++i) {
        infections.push_back (checkpointedInfection (stream));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}

void CommonWithinHost::checkpoint (ostream& stream) {
    WHFalciparum::checkpoint (stream);
    (*pkpdModel) & stream;
    for (std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end(); ++inf) {
        (**inf) & stream;
    }
}
}
}
