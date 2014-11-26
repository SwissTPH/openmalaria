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

#include "Global.h"
#include "WithinHost/DescriptiveWithinHost.h"
#include "WithinHost/Diagnostic.h"
#include "WithinHost/Genotypes.h"
#include "WithinHost/Pathogenesis/PathogenesisModel.h"
#include "util/ModelOptions.h"
#include "PopulationStats.h"
#include "util/StreamValidator.h"
#include "util/errors.h"
#include <cassert>

using namespace std;

namespace OM {
namespace WithinHost {

extern bool bugfix_max_dens;    // DescriptiveInfection.cpp
bool reportPatentInfected = false;

// -----  Initialization  -----

void DescriptiveWithinHostModel::initDescriptive(){
    reportPatentInfected = mon::isUsedM(mon::MHR_PATENT_INFECTIONS);
}

DescriptiveWithinHostModel::DescriptiveWithinHostModel( double comorbidityFactor ) :
        WHFalciparum( comorbidityFactor )
{
    assert( sim::oneTS() == sim::fromDays(5) );
}

DescriptiveWithinHostModel::~DescriptiveWithinHostModel() {}


// -----  Simple infection adders/removers  -----

void DescriptiveWithinHostModel::loadInfection(istream& stream) {
    infections.push_back(DescriptiveInfection(stream));
}

void DescriptiveWithinHostModel::clearInfections( Treatments::Stages stage ){
    for (std::list<DescriptiveInfection>::iterator inf = infections.begin(); inf != infections.end();) {
        if( stage == Treatments::BOTH ||
            (stage == Treatments::LIVER && !inf->bloodStage()) ||
            (stage == Treatments::BLOOD && inf->bloodStage())
        ){
            inf = infections.erase( inf );
        }else{
            ++inf;
        }
    }
    numInfs = infections.size();
}

// -----  Interventions  -----

void DescriptiveWithinHostModel::clearImmunity() {
    for (std::list<DescriptiveInfection>::iterator inf = infections.begin(); inf != infections.end(); ++inf) {
        inf->clearImmunity();
    }
    m_cumulative_h = 0.0;
    m_cumulative_Y_lag = 0.0;
}
void DescriptiveWithinHostModel::importInfection(){
    PopulationStats::totalInfections += 1;
    if( numInfs < MAX_INFECTIONS ){
        PopulationStats::allowedInfections += 1;
        m_cumulative_h += 1;
        numInfs += 1;
        infections.push_back(DescriptiveInfection());
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}


// -----  Density calculations  -----

void DescriptiveWithinHostModel::update(int nNewInfs, vector<double>& genotype_weights,
        double ageInYears, double bsvFactor)
{
    // Cache total density for infectiousness calculations
    assert( Genotypes::N() == 1 );
    m_y_lag.at(sim::ts0().moduloSteps(y_lag_len), 0/*first and only genotype*/) = totalDensity;
    
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    PopulationStats::totalInfections += nNewInfs;
    nNewInfs=min(nNewInfs,MAX_INFECTIONS-numInfs);
    PopulationStats::allowedInfections += nNewInfs;
    numInfs += nNewInfs;
    assert( numInfs>=0 && numInfs<=MAX_INFECTIONS );
    for ( int i=0; i<nNewInfs; ++i ) {
        infections.push_back(DescriptiveInfection());
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
    
    for (std::list<DescriptiveInfection>::iterator inf = infections.begin(); inf != infections.end();) {
        //NOTE: it would be nice to combine this code with that in
        // CommonWithinHost.cpp, but a few changes would be needed:
        // INNATE_MAX_DENS and MAX_DENS_CORRECTION would need to be required
        // (couldn't support old parameterisations using buggy versions of code
        // any more).
        // SP drug action and the PK/PD model would need to be abstracted
        // behind a common interface.
        if ( inf->expired() /* infection has self-terminated */ ||
            (inf->bloodStage() ? treatmentBlood : treatmentLiver) )
        {
            inf=infections.erase(inf);
            numInfs--;
            continue;
        }
        
        // Should be: infStepMaxDens = 0.0, but has some history.
        // See MAX_DENS_CORRECTION in DescriptiveInfection.cpp.
        double infStepMaxDens = timeStepMaxDensity;
        inf->determineDensities(ageInYears, cumulative_h, cumulative_Y, infStepMaxDens, _innateImmSurvFact, bsvFactor);

        if (bugfix_max_dens)
            infStepMaxDens = std::max(infStepMaxDens, timeStepMaxDensity);
        timeStepMaxDensity = infStepMaxDens;

        double density = inf->getDensity();
        totalDensity += density;
        m_cumulative_Y += sim::oneTS().inDays() * density;

        ++inf;
    }
    util::streamValidate( totalDensity );
    assert( (boost::math::isfinite)(totalDensity) );        // inf probably wouldn't be a problem but NaN would be
}


// -----  Summarize  -----

bool DescriptiveWithinHostModel::summarize( const Host::Human& human )const{
    pathogenesisModel->summarize( human );
    
    if( infections.size() > 0 ){
        mon::reportMHI( mon::MHR_INFECTED_HOSTS, human, 1 );
        // (patent) infections are reported by genotype, even though we don't have
        // genotype in this model
        mon::reportMHGI( mon::MHR_INFECTIONS, human, 0, infections.size() );
        if( reportPatentInfected ){
            for (std::list<DescriptiveInfection>::const_iterator inf =
                infections.begin(); inf != infections.end(); ++inf) {
                if( Monitoring::Survey::diagnostic().isPositive( inf->getDensity() ) ){
                    mon::reportMHGI( mon::MHR_PATENT_INFECTIONS, human, 0, 1 );
                }
            }
        }
    }
    
    // Some treatments (simpleTreat with steps=-1) clear infections immediately
    // (and are applied after update()), thus infections.size() may be 0 while
    // totalDensity > 0. Here we report the last calculated density.
    if( Monitoring::Survey::diagnostic().isPositive(totalDensity) ){
        mon::reportMHI( mon::MHR_PATENT_HOSTS, human, 1 );
        mon::reportMHF( mon::MHF_LOG_DENSITY, human, log(totalDensity) );
        return true;    // patent
    }
    return false;       // not patent
}


// -----  Data checkpointing  -----

void DescriptiveWithinHostModel::checkpoint (istream& stream) {
    WHFalciparum::checkpoint (stream);
    for (int i=0; i<numInfs; ++i) {
        loadInfection(stream);  // create infections using a virtual function call
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}
void DescriptiveWithinHostModel::checkpoint (ostream& stream) {
    WHFalciparum::checkpoint (stream);
    BOOST_FOREACH (DescriptiveInfection& inf, infections) {
        inf & stream;
    }
}

char const*const not_impl = "feature not available with the \"descriptive\" within-host model";
void DescriptiveWithinHostModel::treatPkPd(size_t schedule, size_t dosages, double age){
    throw TRACED_EXCEPTION( not_impl, util::Error::WHFeatures ); }

}
}
