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

#include "mon/AgeGroup.h"

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
    // Number of sub-population removals due to first event. Units: humans
    MHR_SUB_POP_REM_FIRST_EVENT,
    // Number of sub-population removals due to expiry. Units: humans
    MHR_SUB_POP_REM_TOO_OLD,
    // Number of infected hosts by genotype. Units: humans
    MHR_INFECTED_GENOTYPE,
    // Number of patent infected hosts by genotype. Units: humans
    MHR_PATENT_GENOTYPE,
    
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
    // Number of pre-erythrocytic vaccine doses deployed. Units: doses (as above)
    MHD_PEV,
    // Number of blood-stage vaccine doses deployed. Units: doses (as above)
    MHD_BSV,
    // Number of transmission-blocking vaccine doses deployed. Units: doses (as above)
    MHD_TBV,
    // Number of bet nets deployed (technically: deployments using the "ITN" model)
    MHD_ITN,
    // Number of IRS spray rounds (technically: deployments using the "IRS" model)
    MHD_IRS,
    // Number of human-vector intervention deployments (technically: deployments using the "GVI" model)
    MHD_GVI,
    // Number of treat intervention deployments (e.g. treatments deployed in an MDA/MSAT campaign)
    MHD_TREAT,
    // Number of screenings done (e.g. tests used in MSAT/T&T)
    MHD_SCREEN,
    // Number of sub-pop recruitments without deployment (deployments of "FIXME" intervention)
    MHD_RECRUIT,
    // Number of deployments (all interventions)
    MHD_ALL_DEPLOYS,
    
    // ———  MHF: measures for human reports (double)  ———
    // Expected number of new infections per human. Units: infections
    MHF_EXPECTED_INFECTED,
    // Report of pyrogenic threshold. Units: ?
    MHF_PYROGENIC_THRESHOLD,
    // Report of log of pyrogenic threshold. Units: ?
    MHF_LOG_PYROGENIC_THRESHOLD,
    // Report of natural log of total parasite density in humans. Units: log(PRBC/μl)
    MHF_LOG_DENSITY,
    // As MHF_LOG_DENSITY, but per genotype
    MHF_LOG_DENSITY_GENOTYPE,
    // Report of age of humans. Units: years
    MHF_AGE,
    
    // ———  MVF: vector (transmission) measures (doubles)  ———
    // Infectiousness of human population to mosquitoes
    MVF_NUM_TRANSMIT,
    // Annual Average Kappa
    MVF_ANN_AVG_K,
    // Input EIR (Expected EIR entered into scenario file). Units: inoculations per adult per time step.
    MVF_INPUT_EIR,
    // Simulated EIR (EIR output by the transmission model). Units: inoculations per adult per time step.
    MVF_SIM_EIR,
    // Total inoculations over survey period per group (age, cohort). Units: inoculations.
    MVF_INOCS,
    // N_v0: emergence of feeding vectors during the last time step. Units: mosquitoes/day
    MVF_LAST_NV0,
    // N_v: vectors seeking to feed during the last time step. Units: mosquitoes/day
    MVF_LAST_NV,
    // O_v: infected vectors seeking to feed during the last time step. Units: mosquitoes/day
    MVF_LAST_OV,
    // S_v: infectious vectors seeking to feed during the last time step. Units: mosquitoes/day
    MVF_LAST_SV,
    
    M_NUM,
    M_OBSOLETE,
    M_ALL_CAUSE_IMR
};

namespace Deploy {
    /// Deployment methods
    enum Method {
        NA = 0,     // not a deployment method
        TIMED = 1<<0,   // mass distribution campaign
        CTS = 1<<1, // continuous deployment (EPI, etc.)
        TREAT = 1<<2    // triggered by case management
    };
}

/// Report some value (integer) to the current survey.
void reportMI( Measure measure, int val );
/// Report some value (floating point) to the current survey.
void reportMF( Measure measure, double val );
/// Report some value (integer) for some human to the current survey.
void reportMHI( Measure measure, const Host::Human& human, int val );
/// Report some value (integer) for some survey, age group and cohort set
void reportMSACI( Measure measure, size_t survey, AgeGroup ageGroup,
                  uint32_t cohortSet, int val );
/// Report some value (integer) for some human and genotype to the current survey.
void reportMHGI( Measure measure, const Host::Human& human, size_t genotype,
                 int val );
/// Report one deployment for some human to the current survey.
void reportMHD( Measure measure, const Host::Human& human,
                Deploy::Method method );

/// Report some value (floating point) for some human to the current survey.
void reportMHF( Measure measure, const Host::Human& human, double val );
/// Report some value (floating point) for the current survey and some age
/// group and cohort set
void reportMACGF( Measure measure, size_t ageIndex, uint32_t cohortSet,
                  size_t genotpye, double val );
/// Report some value (floating point) for some human and genotype to the current survey.
void reportMHGF( Measure measure, const Host::Human& human, size_t genotype,
                 double val );
/// Report some value (floating point) by vector species to the current survey.
void reportMSF( Measure measure, size_t species, double val );
/// Report some value (floating point) by genotype to the current survey.
void reportMSGF( Measure measure, size_t species, size_t genotype, double val );
/// Query whether an output measure is used.
/// This function is not fast, so it is recommended to cache the result.
bool isUsedM( Measure measure );

}
}
#endif
