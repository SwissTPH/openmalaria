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

#ifndef Hmod_CaseManagementCommon
#define Hmod_CaseManagementCommon

// Data that is common to all case-management models

#include "Global.h"
#include "util/AgeGroupInterpolation.h"
#include "Parameters.h"

namespace OM { namespace Clinical {

///@brief Case fatality and sequelae "rate" data
//@{
/** Read community CFR ratio parameter. */
void initCommunityCFR( const Parameters& parameters );

/** Calculate the case fatality "rate" in the community as a function of
 * that in hospitals. */
double getCommunityCFR(double caseFatalityRatio);

/// Age-specific hospital case fatality "rates"
extern util::AgeGroupInterpolator caseFatalityRate;
/// Age-specific in-hospital rates of sequelae given a severe malaria bout
/// Note: out-patients have currently have the same probabilities of sequelae
extern util::AgeGroupInterpolator pSequelaeInpatient;
//@}

} }
#endif
