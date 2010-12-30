/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef H_SurveyMeasure
#define H_SurveyMeasure

/** This file is shared with the BOINC-server apps.
 * It should be kept synchronised somehow...
 *
 * Namespace SM (Survey-Measure or Simulation-Measure) is used to separate
 * contents from standard code base. */
namespace SM {
    
/** Enumeration of reporting options
 *
 * Many are reported per age-group, but to check which actually are you'll have
 * to look through the code.
 * 
 * Don't ever change these names or numbers. The names are used in scenario
 * files, and the numbers in results output/databases. */
enum SurveyMeasure {
    /// Total number of humans
    nHost = 0,
    /// number of infected hosts 
    nInfect = 1,
    /// expected number of infected hosts
    nExpectd= 2,
    /// number of patent hosts
    nPatent= 3,
    /// Sum of the log of the pyrogen threshold
    sumLogPyrogenThres = 4,
    /// Sum of the logarithm of the parasite density
    sumlogDens= 5,
    /// Total infections
    totalInfs= 6,
    /** Infectiousness of human population to mosquitoes
     *
     * Number of hosts transmitting to mosquitoes (i.e. sum of proportion of
     * mosquitoes that get infected). Single value, not per age-group. */
    nTransmit= 7,
    /// Total patent infections
    totalPatentInf= 8,
    /// Contribution to immunity functions
    ///NOTE: not used
    contrib= 9,
    /// Sum of the pyrogenic threshold
    sumPyrogenThresh = 10,
    /// number of treatments (1st line)
    nTreatments1= 11,
    /// number of treatments (2nd line)
    nTreatments2= 12,
    /// number of treatments (inpatient)
    nTreatments3= 13,
    /// number of episodes (uncomplicated)
    nUncomp= 14,
    /// number of episodes (severe)
    nSevere= 15,
    /// cases with sequelae
    nSeq= 16,
    /// deaths in hospital
    nHospitalDeaths= 17,
    /// number of deaths (indirect)
    nIndDeaths= 18,
    /// number of deaths (direct)
    nDirDeaths= 19,
    /// number of EPI vaccine doses given
    nEPIVaccinations= 20,
    /// all cause infant mortality rate
    allCauseIMR= 21,
    /// number of Mass / Campaign vaccine doses given
    nMassVaccinations= 22,
    /// recoveries in hospital
    nHospitalRecovs= 23,
    /// sequelae in hospital
    nHospitalSeqs= 24,
    /// number of IPT Doses
    nIPTDoses= 25,
    /** Annual Average Kappa
     *
     * Calculated once a year as sum of human infectiousness divided by initial
     * EIR summed over a year. Single value, not per age-group. */
    annAvgK= 26,
    /// Number of episodes (non-malaria fever)
    nNMFever= 27,
    /** The total number of inoculations per age group, summed over the
     * reporting period. */
    innoculationsPerAgeGroup = 30,
    
    //BEGIN Per day-of-year data (removed)
    /// Inoculations per human (all ages) per day of year, over the last year.
    /// (Reporting removed.)
    innoculationsPerDayOfYear = 28,
    /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
    /// (Reporting removed.)
    kappaPerDayOfYear = 29,
    //END
    
    /** @brief Vector model parameters.
     *
     * All are point-time outputs, not averages. The Nv0, Nv, Ov and Sv outputs
     * are per-species; the EIR outputs are single values. */
    //@{
    /** Mosquito emergence rate. */
    Vector_Nv0 = 31,
    /// Mosquito population size
    Vector_Nv = 32,
    /// Number of infected mosquitoes
    Vector_Ov = 33,
    /// Number of infectious mosquitoes
    Vector_Sv = 34,
    /** Input EIR (Expected EIR entered into scenario file)
     *
     * Units: inoculations per adult per timestep.
     */
    Vector_EIR_Input = 35,
    /** Simulated EIR (EIR output by the transmission model)
     *
     * Units: inoculations per person per timestep (not per-adult:
     * since children are less available to mosquitoes than adults,
     * this population-average figure ends up being smaller than if
     * all modelled humans were adults).
     */
    Vector_EIR_Simulated = 36,
    //@}
    
    /// @brief EventScheduler reporting (additional to above)
    //@{
    /// Number of Rapid Diagnostic Tests used
    Clinical_RDTs = 39,
    /** Effective total quanty of each drug used orally, in mg.
     * (Per active ingredient abbreviation.)
     * 
     * The quantity is efffective with respect to the cost (see treatment
     * schedule definition). */
    Clinical_DrugUsage = 40,
    /// Direct death on first day of CM (before treatment takes effect)
    Clinical_FirstDayDeaths = 41,
    /// Direct death on first day of CM (before treatment takes effect); hospital only
    Clinical_HospitalFirstDayDeaths = 42,
    //@}
    
    /** The number of actual infections since the last survey. */
    nNewInfections = 43,
    
    /** The number of ITNs delivered by mass distribution since last survey.
     *
     * These are "modelled ITNs": cover only a single person, cannot be passed
     * to someone else for reuse or used for fishing, etc. */
    nMassITNs = 44,
    /** The number of ITNs delivered through EPI since last survey.
     *
     * Comments from nMassITNs apply. */
    nEPI_ITNs = 45,
    /** The number of people newly protected by IRS since last survey.
     *
     * Modelled IRS: affects one person, cannot be plastered over. */
    nMassIRS = 46,
    /** The number of people newly protected by a vector-availability
     * intervention since the last survey. */
    nMassVA = 47,
    
    /// Number of malarial tests via microscopy used
    Clinical_Microscopy = 48,
    /** As Clinical_DrugUsage, but for quatities of drug delivered via IV. */
    Clinical_DrugUsageIV = 49,
    
    /// Number of individuals added to cohort
    nAddedToCohort = 50,
    /// Number of individuals removed from cohort
    nRemovedFromCohort = 51,
    
    /** Number of people (per age group) treated by mass drug administration
     * campaign. (Note that in one day time-step model MDA can be configured
     * as screen-and-treat. This option repeats actual treatments.) */
    nMDAs = 52,
    
    // must be hightest value above plus one
    NUM_SURVEY_OPTIONS	
};
}
#endif
