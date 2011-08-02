/*
  This file is part of OpenMalaria.

  Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

  OpenMalaria is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or (at
  your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "WithinHost/DescriptiveWithinHost.h"
#include "util/ModelOptions.h"
#include "PopulationStats.h"
#include "util/StreamValidator.h"
#include <cassert>

using namespace std;

namespace OM {
namespace WithinHost {

// -----  Initialization  -----

DescriptiveWithinHostModel::DescriptiveWithinHostModel() :
        WithinHostModel()
{
    assert( TimeStep::interval == 5 );
}

DescriptiveWithinHostModel::~DescriptiveWithinHostModel() {
    clearAllInfections();
}


// -----  Simple infection adders/removers  -----

DescriptiveInfection* DescriptiveWithinHostModel::createInfection () {
    return new DescriptiveInfection();
}
void DescriptiveWithinHostModel::loadInfection(istream& stream) {
    infections.push_back(new DescriptiveInfection(stream));
}

void DescriptiveWithinHostModel::clearAllInfections() {
    std::list<DescriptiveInfection*>::iterator inf;
    for (inf=infections.begin(); inf != infections.end(); ++inf) {
        delete *inf;
    }
    infections.clear();
    numInfs=0;
}

// -----  Interventions  -----

void DescriptiveWithinHostModel::immuneSuppression() {
    for (std::list<DescriptiveInfection*>::iterator inf = infections.begin(); inf != infections.end(); ++inf) {
        (*inf)->immuneSuppression();
    }
    _cumulativeh = 0.0;
    _cumulativeYlag = 0.0;
}
void DescriptiveWithinHostModel::importInfection(){
    PopulationStats::totalInfections += 1;
    if( numInfs < MAX_INFECTIONS ){
        PopulationStats::allowedInfections += 1;
        _cumulativeh += 1;
        numInfs += 1;
        infections.push_back(createInfection());
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}


// -----  Density calculations  -----

void DescriptiveWithinHostModel::update(int nNewInfs, double ageInYears, double BSVEfficacy) {
    // Note: adding infections at the beginning of the update instead of the end
    // shouldn't be significant since before latentp delay nothing is updated.
    PopulationStats::totalInfections += nNewInfs;
    nNewInfs=min(nNewInfs,MAX_INFECTIONS-numInfs);
    PopulationStats::allowedInfections += nNewInfs;
    numInfs += nNewInfs;
    assert( numInfs>=0 && numInfs<=MAX_INFECTIONS );
    for ( int i=0; i<nNewInfs; ++i ) {
        infections.push_back(createInfection());
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

    for (std::list<DescriptiveInfection*>::iterator inf = infections.begin(); inf != infections.end();) {
        //NOTE: it would be nice to combine this code with that in
        // CommonWithinHost.cpp, but a few changes would be needed:
        // INNATE_MAX_DENS and MAX_DENS_CORRECTION would need to be required
        // (couldn't support old parameterisations using buggy versions of code
        // any more).
        // IPT model would need to be overhauled; it might be possible to make
        // a pseudo drug model and move SP code there and remove/generalise the
        // rest.
        if ( (*inf)->expired() /* infection too old */
                || eventSPClears(*inf) /* infection cleared by SP in IPT model */
           ) {
            delete *inf;
            inf=infections.erase(inf);
            numInfs--;
            continue;
        }
        
        // Should be: infStepMaxDens = 0.0, but has some history.
        // See MAX_DENS_CORRECTION in DescriptiveInfection.cpp.
        double infStepMaxDens = timeStepMaxDensity;
        (*inf)->determineDensities(ageInYears, cumulativeh, cumulativeY, infStepMaxDens, _innateImmSurvFact, BSVEfficacy);

        IPTattenuateAsexualDensity (*inf);

        if (util::ModelOptions::option (util::MAX_DENS_CORRECTION))
            infStepMaxDens = std::max(infStepMaxDens, timeStepMaxDensity);
        timeStepMaxDensity = infStepMaxDens;

        totalDensity += (*inf)->getDensity();
        (*inf)->determineDensityFinal ();
        _cumulativeY += TimeStep::interval*(*inf)->getDensity();

        ++inf;
    }
    assert( numInfs == static_cast<int>(infections.size()) );
    util::streamValidate( totalDensity );
    assert( totalDensity == totalDensity );        // inf probably wouldn't be a problem but NaN would be

    IPTattenuateAsexualMinTotalDensity();
}


// -----  Summarize  -----

int DescriptiveWithinHostModel::countInfections (int& patentInfections) {
    if (infections.empty()) return 0;
    patentInfections = 0;
    for (std::list<DescriptiveInfection*>::iterator inf=infections.begin();
            inf != infections.end(); ++inf) {
        if ((*inf)->getDensity() > detectionLimit)
            patentInfections++;
    }
    return infections.size();
}


// -----  Data checkpointing  -----

void DescriptiveWithinHostModel::checkpoint (istream& stream) {
    WithinHostModel::checkpoint (stream);
    for (int i=0; i<numInfs; ++i) {
        loadInfection(stream);  // create infections using a virtual function call
    }
}
void DescriptiveWithinHostModel::checkpoint (ostream& stream) {
    WithinHostModel::checkpoint (stream);
    BOOST_FOREACH (DescriptiveInfection* inf, infections) {
        (*inf) & stream;
    }
}

}
}
