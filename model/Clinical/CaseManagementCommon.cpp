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

#include "Clinical/CaseManagementCommon.h"
#include "util/checkpoint_containers.h"

namespace OM { namespace Clinical {

//log odds ratio of case-fatality in community compared to hospital
double oddsRatioThreshold;

util::AgeGroupInterpolator caseFatalityRate;
util::AgeGroupInterpolator pSequelaeInpatient;

///@brief Infant death summaries (checkpointed).
//@{
vector<int> infantDeaths;
vector<int> infantIntervalsAtRisk;
//@}

/// Non-malaria mortality in under 1year olds.
/// Set by init ()
double nonMalariaMortality;


void initCMCommon( const OM::Parameters& parameters ){
    oddsRatioThreshold = exp( parameters[Parameters::LOG_ODDS_RATIO_CF_COMMUNITY] );
    infantDeaths.resize(TimeStep::stepsPerYear);
    infantIntervalsAtRisk.resize(TimeStep::stepsPerYear);
    nonMalariaMortality=parameters[Parameters::NON_MALARIA_INFANT_MORTALITY];
}

void mainSimInitCMCommon () {
    for (TimeStep i(0);i<TimeStep::intervalsPerYear; ++i) {
        Clinical::infantIntervalsAtRisk[i.asInt()]=0;
        Clinical::infantDeaths[i.asInt()]=0;
    }
}

void staticCheckpointCMCommon (istream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
}
void staticCheckpointCMCommon (ostream& stream) {
    infantDeaths & stream;
    infantIntervalsAtRisk & stream;
}


double getCommunityCFR (double caseFatalityRatio){
    double x = caseFatalityRatio * oddsRatioThreshold;
    return x / (1 - caseFatalityRatio + x);
}



double infantAllCauseMort(){
    double infantPropSurviving=1.0;       // use to calculate proportion surviving
    for (TimeStep i(0);i<TimeStep::intervalsPerYear; ++i) {
        // multiply by proportion of infants surviving at each interval
        infantPropSurviving *= double(infantIntervalsAtRisk[i.asInt()] - infantDeaths[i.asInt()])
            / double(infantIntervalsAtRisk[i.asInt()]);
    }
    // Child deaths due to malaria (per 1000), plus non-malaria child deaths. Deaths per 1000 births is the return unit.
    return (1.0 - infantPropSurviving) * 1000.0 + nonMalariaMortality;
}

} }
