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

#ifndef Hmod_CaseManagementCommon
#define Hmod_CaseManagementCommon

// Data that is common to all case-management models

#include "Global.h"
#include "util/AgeGroupInterpolation.h"
#include "Parameters.h"
#include <vector>

namespace OM { namespace Clinical {

/** Initialise parameters. */
void initCMCommon( const Parameters& parameters, SimTime hsMemory );

void mainSimInitCMCommon ();

/// Static checkpointing
void staticCheckpointCMCommon (istream& stream);
void staticCheckpointCMCommon (ostream& stream);

/// True if bug-fix is enabled; do not set externally
extern bool indirectMortBugfix;

/** The maximum age of a sickness bout, for
 * another bout to be considered part of the same episode.
 * 
 * Used by both the clinical models in roughly the same way, but will have
 * different values in each to match Global::interval. */
extern SimTime healthSystemMemory;


///@brief Case fatality and sequelae "rate" data
//@{
/** Calculate the case fatality "rate" in the community as a function of
 * that in hospitals. */
double getCommunityCFR(double caseFatalityRatio);

/// Age-specific hospital case fatality "rates"
extern util::AgeGroupInterpolator caseFatalityRate;
/// Age-specific in-hospital rates of sequelae given a severe malaria bout
/// Note: out-patients have currently have the same probabilities of sequelae
extern util::AgeGroupInterpolator pSequelaeInpatient;
//@}

/** Calculate infant mortality as deaths/1000 livebirths for the whole main-
 * simulation period (not as deaths/1000 years-at-risk per survey).
 * 
 * This mimicks field data on all-cause mortality in infants.
 * Uses the kaplan-meier method because the demography was set up to provide
 * a stable age-distribution but unfortunately does not accurately describe
 * death rates. The kaplan-meier estimate is the product of the proportion of
 * infants survivng at each interval. */
double infantAllCauseMort();

///@brief Statistics, updated directly by ClinicalModel::updateInfantDeaths
//@{
extern std::vector<int> infantDeaths;
extern std::vector<int> infantIntervalsAtRisk;
//@}

} }
#endif
