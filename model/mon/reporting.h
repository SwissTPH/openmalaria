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

#ifndef H_OM_mon_reporting
#define H_OM_mon_reporting

#include "Monitoring/AgeGroup.h"
#include <boost/integer_traits.hpp>

/** This header handles reporting of data and querying of which outputs are
 * active.
 *
 * It does not handle reading configuration or writing to output files. */
namespace OM {
namespace Host {
    class Human;
}
namespace mon {

/** This enum lists monitoring measures.
 *
 * It does not directly correspond to output codes but rather to things that
 * the model can report.
*/
enum Measure{
    // ———  MHR: measures for human reports (integers)  ———
    // Number of hosts. Units: humans
    MHR_HOSTS,
    // Number of infected hosts. Units: humans
    MHR_INFECTED_HOSTS,
    // Number of patent infected hosts. Units: humans
    MHR_PATENT_HOSTS,
    // Number of infections in humans. Units: infections
    MHR_INFECTIONS,
    // Number of patent infections in humans. Units: infections
    MHR_PATENT_INFECTIONS,
    // Number of new infections in humans. Units: infections
    MHR_NEW_INFECTIONS,
    
    // ———  MHT: measures for human treatments (integers)  ———
    // Number of first line treatments in humans. Units: treatments (whole courses)
    MHT_TREATMENTS_1,
    // Number of second line treatments in humans. Units: treatments (whole courses)
    MHT_TREATMENTS_2,
    // Number of severe/in-hospital treatments in humans. Units: treatments (whole courses)
    MHT_TREATMENTS_3,
    // Number of treatments for non-malaria infections. Units: treatments (whole courses)
    MHT_NMF_TREATMENTS /* also known as antibiotics */,
    // Number of treatments of primaquine. Units: treatments (whole courses)
    MHT_PQ_TREATMENTS,
    // Number of diagnostics used during treatment. Units: diagnostics
    MHT_TREAT_DIAGNOSTICS,
    
    // ———  MHE: measures for human episodes (integers)  ———
    // Number of uncomplicated fever episodes in humans. Units: cases
    MHE_UNCOMPLICATED_EPISODES,
    // Number of severe fever episodes in humans. Units: cases
    MHE_SEVERE_EPISODES,
    // Number of fever episodes in humans not due to malaria. Units: cases
    MHE_NON_MALARIA_FEVERS,
    
    // ———  MHO: outcomes  ———
    // Number of human patients dying directly due to malaria. Units: cases
    MHO_DIRECT_DEATHS,
    // Number of human patients dying indirectly (delayed deaths) due to malaria. Units: cases
    MHO_INDIRECT_DEATHS,
    // Number of human patients recovering with sequelae. Units: cases
    MHO_SEQUELAE,
    // Number of human patients dying in hospital (directly) due to malaria. Units: cases
    MHO_HOSPITAL_DEATHS,
    // Number of human patients fully recovering in hospital. Units: cases
    MHO_HOSPITAL_RECOVERIES,
    // Number of human patients recovering with sequelae in hospital. Units: cases
    MHO_HOSPITAL_SEQUELAE,
    // Number of human patients dying as a direct result of non-malaria fever. Units: cases
    MHO_NMF_DEATHS,
    // Number of human patients dying on the first day of the episode due to malaria. Units: cases
    MHO_FIRST_DAY_DEATHS,
    // Number of human patients dying on their first day in hospital due to malaria. Units: cases
    MHO_HOSPITAL_FIRST_DAY_DEATHS,
    
    // ———  MHD: intervention deployments  ———
    /** Number of vaccine doses deployed. Units: doses (including first dose,
     * second dose, booster doses, etc.).
     * 
     * Since schema 22, each vaccine type may be deployed independently. To be
     * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
     * (with an "effect" element) will be reported. */
    MHD_VACCINATIONS,
    //TODO: also more specific outputs?
    
    // ———  MHF: measures for human reports (double)  ———
    // Expected number of new infections per human. Units: infections
    MHF_EXPECTED_INFECTED,
    // Report of log of pyrogenic threshold. Units: ?
    MHF_LOG_PYROGENIC_THRESHOLD,
    // Report of natural log of total parasite density in humans. Units: log(PRBC/μl)
    MHF_LOG_DENSITY,
    // Report of pyrogenic threshold. Units: ?
    MHF_PYROGENIC_THRESHOLD,
    // Report of age of humans. Units: years
    MHF_AGE,
    
    M_NUM
};

/// For surveys and measures to say something shouldn't be reported
const size_t NOT_USED = boost::integer_traits<size_t>::const_max;

/// The current survey number (read only usage)
extern size_t currentSurvey;

namespace Deploy {
    /// Deployment methods
    enum Method {
        NA = 0,     // not a deployment method
        TIMED = 1<<0,   // mass distribution campaign
        CTS = 1<<1, // continuous deployment (EPI, etc.)
        TREAT = 1<<2    // triggered by case management
    };
}

/// Report some value (integer) for some human to the current survey.
void reportMHI( Measure measure, const Host::Human& human, int val );
/// Report some value (floating point) for some human to the current survey.
void reportMHF( Measure measure, const Host::Human& human, double val );
/// Report some value (integer) for some survey, age group and cohort set
void reportMSACI( Measure measure, size_t survey, Monitoring::AgeGroup ageGroup, uint32_t cohortSet, int val );
/// Report one deployment for some human to the current survey.
void reportMHD( Measure measure, const Host::Human& human, Deploy::Method method );
/// Query whether an output measure is used.
bool isUsedM( Measure measure );

}
}
#endif
