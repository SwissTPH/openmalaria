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

#include <cassert>
#include "WithinHost/DescriptiveWithinHost.h"
#include "util/ModelOptions.h"
#include "PopulationStats.h"
#include "util/StreamValidator.h"

using namespace std;

namespace OM { namespace WithinHost {

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

void DescriptiveWithinHostModel::newInfection(TimeStep now) {
    ++PopulationStats::totalInfections;
    if (numInfs < MAX_INFECTIONS) {
        infections.push_back(new DescriptiveInfection(now));
        numInfs++;
        ++PopulationStats::allowedInfections;
    }
    assert( numInfs == static_cast<int>(infections.size()) );
}
void DescriptiveWithinHostModel::loadInfection(istream& stream) {
    infections.push_back(new DescriptiveInfection(stream));
}

void DescriptiveWithinHostModel::clearAllInfections() {
    std::list<DescriptiveInfection*>::iterator i;
    for (i=infections.begin(); i != infections.end(); ++i) {
        delete *i;
    }
    infections.clear();
    numInfs=0;
}

// -----  Interventions  -----
void DescriptiveWithinHostModel::immuneSuppression() {
    for (std::list<DescriptiveInfection*>::iterator it = infections.begin(); it != infections.end(); ++it) {
        (*it)->immuneSuppression();
    }
    _cumulativeh = 0.0;
    _cumulativeYlag = 0.0;
}


// -----  Density calculations  -----

void DescriptiveWithinHostModel::calculateDensities(double ageInYears, double BSVEfficacy) {
    updateImmuneStatus ();

    totalDensity = 0.0;
    timeStepMaxDensity = 0.0;

    // Values of _cumulativeh/Y at beginning of step
    // (values are adjusted for each infection)
    double cumulativeh=_cumulativeh;
    double cumulativeY=_cumulativeY;

    std::list<DescriptiveInfection*>::iterator iter=infections.begin();
    while (iter != infections.end()) {
        if ( (*iter)->expired() /* infection too old */
                || eventSPClears(*iter) /* infection cleared by SP in IPT model */
           ) {
            delete *iter;
            iter=infections.erase(iter);
            numInfs--;
        }
        else {
            // Should be: infStepMaxDens = 0.0, but has some history.
            // See MAX_DENS_CORRECTION in DescriptiveInfection.cpp.
            double infStepMaxDens = timeStepMaxDensity;
            (*iter)->determineDensities(ageInYears, cumulativeh, cumulativeY, infStepMaxDens, _innateImmSurvFact, BSVEfficacy);

            IPTattenuateAsexualDensity (*iter);

            if (util::ModelOptions::option (util::MAX_DENS_CORRECTION))
                infStepMaxDens = std::max(infStepMaxDens, timeStepMaxDensity);
            timeStepMaxDensity = infStepMaxDens;

            totalDensity += (*iter)->getDensity();
            if ((*iter)->getStartDate() == TimeStep::simulation1()-TimeStep(1)) {
                _cumulativeh++;
            }
            (*iter)->determineDensityFinal ();
            _cumulativeY += TimeStep::interval*(*iter)->getDensity();

            ++iter;
        }
    }
    assert( numInfs == static_cast<int>(infections.size()) );

    IPTattenuateAsexualMinTotalDensity();
    util::streamValidate( totalDensity );
}


// -----  Summarize  -----

int DescriptiveWithinHostModel::countInfections (int& patentInfections) {
    if (infections.empty()) return 0;
    patentInfections = 0;
    for (std::list<DescriptiveInfection*>::iterator iter=infections.begin();
            iter != infections.end(); ++iter) {
        if ((*iter)->getDensity() > detectionLimit)
            patentInfections++;
    }
    return infections.size();
}


// -----  Data checkpointing  -----

void DescriptiveWithinHostModel::checkpoint (istream& stream) {
    WithinHostModel::checkpoint (stream);
    for (int i=0; i<numInfs; ++i) {
        loadInfection(stream);	// create infections using a virtual function call
    }
}
void DescriptiveWithinHostModel::checkpoint (ostream& stream) {
    WithinHostModel::checkpoint (stream);
    BOOST_FOREACH (DescriptiveInfection* inf, infections) {
        (*inf) & stream;
    }
}

} }
