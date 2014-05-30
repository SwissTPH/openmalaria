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

using namespace std;

namespace OM {
namespace WithinHost {

CommonInfection* (* CommonWithinHost::createInfection) (uint32_t protID);
CommonInfection* (* CommonWithinHost::checkpointedInfection) (istream& stream);


// -----  Initialization  -----

CommonWithinHost::CommonWithinHost( double comorbidityFactor ) :
        WHFalciparum( comorbidityFactor ), pkpdModel(PkPd::PkPdModel::createPkPdModel ())
{
    assert( TimeStep::interval == 1 || TimeStep::interval == 5 );
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

void CommonWithinHost::medicate(string drugName, double qty, double time, double duration, double bodyMass) {
    pkpdModel->medicate(drugName, qty, time, duration, bodyMass);
}
void CommonWithinHost::clearImmunity() {
    for (std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end(); ++inf) {
        (*inf)->clearImmunity();
    }
    _cumulativeh = 0.0;
    _cumulativeYlag = 0.0;
}
void CommonWithinHost::importInfection(){
    PopulationStats::totalInfections += 1;
    if( numInfs < MAX_INFECTIONS ){
        PopulationStats::allowedInfections += 1;
        _cumulativeh += 1;
        numInfs += 1;
        infections.push_back(createInfection (pkpdModel->new_proteome_ID ()));
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}


// -----  Density calculations  -----

void CommonWithinHost::update(int nNewInfs, double ageInYears, double bsvFactor) {
    // Cache total density for infectiousness calculations
    _ylag[mod_nn(TimeStep::simulation.asInt(),_ylagLen)] = totalDensity;
    
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
    
    // As in AJTMH p22, cumulativeh (X_h + 1) doesn't include infections added
    // this time-step and cumulativeY only includes past densities.
    double cumulativeh=_cumulativeh;
    double cumulativeY=_cumulativeY;
    _cumulativeh += nNewInfs;

    bool treatmentLiver = treatExpiryLiver >= TimeStep::simulation;
    bool treatmentBlood = treatExpiryBlood >= TimeStep::simulation;
    bool treatmentBoth = treatmentLiver && treatmentBlood;
    
    for (std::list<CommonInfection*>::iterator inf = infections.begin(); inf != infections.end();) {
        double survivalFactor = bsvFactor * _innateImmSurvFact;
        survivalFactor *= (*inf)->immunitySurvivalFactor(ageInYears, cumulativeh, cumulativeY);
        survivalFactor *= pkpdModel->getDrugFactor((*inf)->get_proteome_ID());
        
        for( int step = 0, steps = TimeStep::interval; step < steps; ++step ){
            // Note: this is only one treatment model; there is also the PK/PD model
            bool expires = treatmentBoth /* effective treatment targeting both stages */ ||
                (treatmentBlood && (*inf)->bloodStage()) /* blood stage treatment */ ||
                (treatmentLiver && !(*inf)->bloodStage()) /* liver stage treatment */;
            
            if( !expires ) /* no expiry due to simple treatment model; do update */
                expires = (*inf)->update(survivalFactor) /* returns true on expiry (due to PK/PD or self-termination) */;
            
            if( expires ){
                delete *inf;
                inf = infections.erase(inf);        // inf points to next infection now so don't increment with ++inf
                --numInfs;
                goto CONTINUE_OUTER;   // infection no longer exists so skip the rest
            }
            
            double density = (*inf)->getDensity();
            totalDensity += density;
            timeStepMaxDensity = max(timeStepMaxDensity, density);
            _cumulativeY += density;
        }
        
        ++inf;
        CONTINUE_OUTER:; // yes, it's a hideous goto — but C++98 doesn't have another way of continuing the outer loop
    }
    util::streamValidate(totalDensity);
    assert( (boost::math::isfinite)(totalDensity) );        // inf probably wouldn't be a problem but NaN would be

    pkpdModel->decayDrugs ();
}

void CommonWithinHost::addProphylacticEffects(const vector<double>& pClearanceByTime) {
    // this should actually be easy; it just isn't needed yet
    throw util::unimplemented_exception( "prophylactic effects on 1-day timestep" );
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
