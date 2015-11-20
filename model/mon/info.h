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

#ifndef H_OM_mon_info
#define H_OM_mon_info

#include "Global.h"     // SimTime

#include <boost/integer_traits.hpp>

/** This header provides information from the reporting system. */
namespace OM {
namespace mon {
// Not 'private' but still not for use externally:
namespace impl {
    // Consts (set during program start-up):
    extern size_t nSurveys;     // number of reported surveys
    extern size_t nCohortSets;
    // Variables (checkpointed):
    extern bool isInit; // set true after "initialisation" survey at intervention time 0
    extern size_t survNumEvent, survNumStat;
    extern SimTime nextSurveyTime;
}

/// For surveys and measures to say something shouldn't be reported
const size_t NOT_USED = boost::integer_traits<size_t>::const_max;

/** Line end character. Use Unix line endings to save a little size. */
const char lineEnd = '\n';

/// The current survey number (can be passed back to 'event' report functions taking
/// survey times). May have the special value NOT_USED.
inline size_t eventSurveyNumber(){ return impl::survNumEvent; }

/// Whether the current survey is reported.
/// 
/// Exception: there is a dummy survey at intervention time 0 which is not
/// reported but acts like it is to set survey variables.
inline bool isReported(){ return !impl::isInit || impl::survNumStat != NOT_USED; }

/** Time the current (next) survey ends at, or sim::never() if no more
 * surveys take place. */
SimTime nextSurveyTime();

/** Return the time of the final survey.
 *
 * We use this to control when the simulation ends.
 * This isn't quite the same as before when the simulation end was
 * explicitly specified and has a small affect on
 * infantAllCauseMortality (survey 21) output. */
SimTime finalSurveyTime();

/// The number of cohort sets
inline size_t numCohortSets(){ return impl::nCohortSets; }

/** Humans should store a "cohort set" identifier which is initially 0.
 * Whenever a human gains or loses membership status in some
 * sup-population, it should update that value with this function.
 * 
 * @param old       Old identifier value (initially 0)
 * @param subPop    Sub-population to which membership status changed
 * @param isMember  New membership status
 * @returns         New identifier value */
uint32_t updateCohortSet( uint32_t old, interventions::ComponentId subPop,
        bool isMember );

} }
#endif
