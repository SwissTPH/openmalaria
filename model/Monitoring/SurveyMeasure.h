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
    BLANK = 0,  // temporary: items moved to OutputMeasures.hpp
    
    /** Infectiousness of human population to mosquitoes
     *
     * Number of hosts transmitting to mosquitoes (i.e. proportion of
     * mosquitoes that get infected multiplied by human population size).
     * Single value, not per age-group. */
    nTransmit= 7,
//     /// Contribuion to immunity functions (removed)
//     contrib= 9,
    
    /** Number of vaccine doses given via EPI.
     * 
     * Since schema 22, each vaccine type may be deployed independently. To be
     * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
     * (with an "effect" element) will be reported. */
    nEPIVaccinations= 20,
    
    /** All cause infant mortality rate
     * 
     * Reports death rate of infants due to all causes (malaria as modelled
     * plus fixed non-malaria attribution). Calculated via Kaplan-Meier method.
     * Units: deaths per thousand births.
     * 
     * For more info, see ClinicalModel::infantAllCauseMort() (ClinicalModel.h,
     * line 74).
     */
    allCauseIMR= 21,
    
    /** Number of vaccine doses given via mass campaign.
     * 
     * Since schema 22, each vaccine type may be deployed independently. To be
     * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
     * (with an "effect" element) will be reported. */
    nMassVaccinations= 22,
    
//     /// Number of IPT Doses (removed together with IPT model)
//     nIPTDoses= 25,
    
    /** Annual Average Kappa
     *
     * Calculated once a year as sum of human infectiousness divided by initial
     * EIR summed over a year. Single value, not per age-group. */
    annAvgK= 26,
    
    /** The total number of inoculations per age group, summed over the
     * reporting period. */
    innoculationsPerAgeGroup = 30,
    
    //BEGIN Per day-of-year data (removed)
//     /// Inoculations per human (all ages) per day of year, over the last year.
//     /// (Reporting removed.)
//     innoculationsPerDayOfYear = 28,
//     /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
//     /// (Reporting removed.)
//     kappaPerDayOfYear = 29,
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
     * Units: inoculations per adult per time step.
     */
    inputEIR = 35,
    /** Simulated EIR (EIR output by the transmission model)
     *
     * Units: inoculations per person per time step (not per-adult:
     * since children are less available to mosquitoes than adults,
     * this population-average figure ends up being smaller than if
     * all modelled humans were adults).
     */
    simulatedEIR = 36,
    //@}
    
    /// @brief EventScheduler reporting (additional to above)
    //@{
    /// Number of Rapid Diagnostic Tests used
    Clinical_RDTs = 39,
    /* Effective total quanty of each drug used orally, in mg.
     * (Per active ingredient abbreviation.)
     * 
     * The quantity is efffective with respect to the cost (see treatment
     * schedule definition).
     * 
     * Reporting removed. */
    //Clinical_DrugUsage = 40,
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
    /** Defunct; was used by "vector availability" intervention (which is now a
     * sub-set of GVI). */
//     nMassVA = 47,
    
    /// Number of malarial tests via microscopy used
    Clinical_Microscopy = 48,
    /* As Clinical_DrugUsage, but for quatities of drug delivered via IV. */
    //Clinical_DrugUsageIV = 49,
    
//     /// Number of cohort recruitments removed)
//     nAddedToCohort = 50,
//     /// Number of individuals removed from cohort (removed)
//     nRemovedFromCohort = 51,
    
    /** Number of people (per age group) treated by mass drug administration
     * campaign. (Note that in one day time-step model MDA can be configured
     * as screen-and-treat. This option reports treatments administered not
     * the number of tests used.) */
    nMDAs = 52,
    
//     /// Number of antibiotic treatments given (disabled — not used)
//     nAntibioticTreatments = 54,
    
    /** Report the number of screenings used in a mass screen-and-treat
     * operation. */
    nMassScreenings = 55,
    
    /** Report the number of mass deployments of generic vector interventions.
     * 
     * Note: this is a provisionary reporting measure. Like many other measures,
     * it is insufficient now that multiple descriptions of any intervention
     * type are possible. */
    nMassGVI = 56,
    
    /** Number of IRS deployments via continuous deployment. */
    nCtsIRS = 57,
    
    /** Number of GVI deployments via continuous deployment. */
    nCtsGVI = 58,
    
    /** Number of "MDA" deployments via continuous deployment.
     * 
     * Note: MDA stands for mass drug administration, but the term has come to
     * be used more flexibly by OpenMalaria, including optional screening and
     * deployment through age-based systems. */
    nCtsMDA = 59,
    
    /** Number of diagnostics used by "MDA" distribution through continuous
     * methods. Can be higher than nCtsMDA since drugs are administered only
     * when the diagnostic is positive. Also see nCtsMDA description. */
    nCtsScreenings = 60,
    
    /** Number of removals from a sub-population due to expiry of duration of
     * membership (e.g. intervention too old). */
    nSubPopRemovalTooOld = 61,
    /** Number of removals from a sub-population due to first
     * infection/bout/treatment (see onFirstBout & co). */
    nSubPopRemovalFirstEvent = 62,
    
    /** Report the number of Primaquine treatments given. */
    nPQTreatments = 63,
    
    /** Report the number of diagnostics used during treatment.
     * 
     * This is not the same as Clinical_RDTs + Clinical_Microscopy: those
     * outputs are used by the "event scheduler" 1-day time step clinical
     * model, whereas this output is used by the 5-day time step model. */
    nTreatDiagnostics = 64,
    
    /** Number of "recruitment only" recruitments via timed deployment. */
    nMassRecruitOnly = 65,
    /** Number of "recruitment only" recruitments via age-based deployment. */
    nCtsRecruitOnly = 66,
    
    /** Number of deployments (of all intervention components) triggered by
     * treatment (case management). */
    nTreatDeployments = 67,
    
    // must be hightest value above plus one
    NUM_SURVEY_OPTIONS	
};
}
#endif
