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

#ifndef H_OM_mon_cpp
// This is the only way we need to use this header:
#error "Only include from mon.cpp"
#endif

#include "mon/reporting.h"      // for Measure enum
#include <map>
#include <string>

namespace OM {
namespace mon {

// Describes each "measure" to be output
struct OutMeasure{
    // Number used in output to identify this measure/aggregation.
    // This is what identifies this "measure".
    int outId;
    
    // The following control what is reported by this measure:
    // Set m to M_NUM for obsolete/special outputs
    Measure m;  // internal measure (e.g. MHR_HOSTS) this comes from
    bool isDouble;  // false: type is int; true: type is double
    bool byAge; // segregate by age
    bool byCohort;      // segregate by cohort
    bool bySpecies;     // segregate by species of vector
    bool byGenotype;    // segregate by genotype of parasite
    uint8_t method;     // deployment method (see above)
    
    // Convenience constructors:
    OutMeasure() : outId(-1), m(M_NUM), isDouble(false), byAge(false),
                byCohort(false), bySpecies(false), byGenotype(false), method(0) {}
    OutMeasure( int outId, Measure m, bool isDouble, bool byAge, bool byCohort,
                bool bySpecies, bool byGenotype, uint8_t method ) :
        outId(outId), m(m), isDouble(isDouble), byAge(byAge), byCohort(byCohort),
        bySpecies(bySpecies), byGenotype(byGenotype), method(method) {}
    // Simple reports
    static OutMeasure value( int outId, Measure m, bool isDouble ){
        return OutMeasure( outId, m, isDouble, false, false, false, false, Deploy::NA );
    }
    // Something with reports segregated by human age and cohort membership
    static OutMeasure humanAC( int outId, Measure m, bool isDouble ){
        return OutMeasure( outId, m, isDouble, true, true, false, false, Deploy::NA );
    }
    // Something with reports segregated by human age, cohort membership
    // and parasite genotype.
    static OutMeasure humanACG( int outId, Measure m, bool isDouble ){
        return OutMeasure( outId, m, isDouble, true, true, false, true, Deploy::NA );
    }
    // Reports by mosquito species and optionally parasite genotype.
    // All are floating point (currently).
    static OutMeasure species( int outId, Measure m, bool byGenotype ){
        return OutMeasure( outId, m, true, false, false, true, byGenotype, Deploy::NA );
    }
    // Deployments with reports segregated by human age and cohort membership
    // Method can be Deploy::NA to not match deployments (but in this case,
    // better to use a different constructor), or it can be one of the three
    // deployment methods (to only count reports of that type of deployment),
    // or it can be a bit-or-ed combination of any of the three methods (to
    // count deployments of multiple types simultaneously).
    static OutMeasure humanDeploy( int outId, Measure m, Deploy::Method method ){
        assert( method >= 0 && method <= (Deploy::TIMED|Deploy::CTS|Deploy::TREAT) );
        return OutMeasure( outId, m, false, true, true, false, false, method );
    }
    static OutMeasure obsolete( int outId ){
        return OutMeasure( outId, M_OBSOLETE, false, false, false, false, false, Deploy::NA );
    }
};

// These are all output measures set by a name in the XML
// Example: nHosts
typedef std::map<std::string,OutMeasure> NamedMeasureMapT;
NamedMeasureMapT namedOutMeasures;

// This method defines output measures accepted by name in the XML (e.g.
// "nHost") and their numeric output identifier (i.e. measure column of
// outputs), type of output (integer or floating point), aggregation, and the
// corresponding internal measure code.
void defineOutMeasures(){
    //NOTE: measures are ordered by their output codes.
    // Add new outputs with next available code at end of list.
    // Don't reuse old codes.
    
    /// Total number of humans
    namedOutMeasures["nHost"] = OutMeasure::humanAC( 0, MHR_HOSTS, false );
    /** The number of human hosts with an infection (patent or not) at the time
     * the survey is taken. */
    namedOutMeasures["nInfect"] = OutMeasure::humanAC( 1, MHR_INFECTED_HOSTS, false );
    /** Expected number of infected hosts
     * 
     * This is the sum of the probabilities, across all time steps since the
     * last survey, of each host becoming infected on that time step. */
    namedOutMeasures["nExpectd"] =
        OutMeasure::humanAC( 2, MHF_EXPECTED_INFECTED, true );
    /** The number of human hosts whose total (blood-stage) parasite density is
     * above the detection threshold */
    namedOutMeasures["nPatent"] = OutMeasure::humanAC( 3, MHR_PATENT_HOSTS, false );
    /// Sum of log(1 + p) where p is the pyrogenic threshold
    namedOutMeasures["sumLogPyrogenThres"] =
        OutMeasure::humanAC( 4, MHF_LOG_PYROGENIC_THRESHOLD, true );
    /** Sum (across hosts) of the natural logarithm of the parasite density of
     * hosts with detectable parasite density (patent according to the
     * monitoring diagnostic). */
    namedOutMeasures["sumlogDens"] = OutMeasure::humanAC( 5, MHF_LOG_DENSITY, true );
    /** The total number of infections in the population: includes both blood
     * and liver stages. Vivax: this is the number of broods. */
    namedOutMeasures["totalInfs"] = OutMeasure::humanACG( 6, MHR_INFECTIONS, false );
    /** Infectiousness of human population to mosquitoes
     *
     * Number of hosts transmitting to mosquitoes (i.e. proportion of
     * mosquitoes that get infected multiplied by human population size).
     * Single value, not per age-group. */
    namedOutMeasures["nTransmit"] = OutMeasure::value( 7, MVF_NUM_TRANSMIT, true );
    /** The sum of all detectable infections (where blood stage parasite
     * density is above the detection limit) across all human hosts.
     * Vivax: the number of broods with an active blood stage. */
    namedOutMeasures["totalPatentInf"] =
        OutMeasure::humanACG( 8, MHR_PATENT_INFECTIONS, false );
    /// Contribuion to immunity functions (removed)
    namedOutMeasures["contrib"] = OutMeasure::obsolete( 9 );
    /// Sum of the pyrogenic threshold
    namedOutMeasures["sumPyrogenThresh"] =
        OutMeasure::humanAC( 10, MHF_PYROGENIC_THRESHOLD, true );
    /// number of treatments (1st line)
    namedOutMeasures["nTreatments1"] = OutMeasure::humanAC( 11, MHT_TREATMENTS_1, false );
    /// number of treatments (2nd line)
    namedOutMeasures["nTreatments2"] = OutMeasure::humanAC( 12, MHT_TREATMENTS_2, false );
    /// number of treatments (inpatient)
    namedOutMeasures["nTreatments3"] = OutMeasure::humanAC( 13, MHT_TREATMENTS_3, false );
    /// number of episodes (uncomplicated)
    namedOutMeasures["nUncomp"] =
        OutMeasure::humanAC( 14, MHE_UNCOMPLICATED_EPISODES, false );
    /// number of episodes (severe)
    namedOutMeasures["nSevere"] =
        OutMeasure::humanAC( 15, MHE_SEVERE_EPISODES, false );
    /// cases with sequelae
    namedOutMeasures["nSeq"] = OutMeasure::humanAC( 16, MHO_SEQUELAE, false );
    /// deaths in hospital
    namedOutMeasures["nHospitalDeaths"] =
        OutMeasure::humanAC( 17, MHO_HOSPITAL_DEATHS, false );
    /// Number of deaths indirectly caused by malaria
    namedOutMeasures["nIndDeaths"] =
        OutMeasure::humanAC( 18, MHO_INDIRECT_DEATHS, false );
    /// Number of deaths directly caused by malaria
    namedOutMeasures["nDirDeaths"] =
        OutMeasure::humanAC( 19, MHO_DIRECT_DEATHS, false );
    /** Number of vaccine doses given via EPI.
     * 
     * Since schema 22, each vaccine type may be deployed independently. To be
     * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
     * (with an "effect" element) will be reported. */
    namedOutMeasures["nEPIVaccinations"] =
        OutMeasure::humanDeploy( 20, MHD_VACCINATIONS, Deploy::CTS );
    /** All cause infant mortality rate
     * 
     * Reports death rate of infants due to all causes (malaria as modelled
     * plus fixed non-malaria attribution). Calculated via Kaplan-Meier method.
     * Units: deaths per thousand births. */
    namedOutMeasures["allCauseIMR"] =
        OutMeasure::value( 21, M_ALL_CAUSE_IMR, true );
    /** Number of vaccine doses given via mass campaign.
     * 
     * Since schema 22, each vaccine type may be deployed independently. To be
     * roughly backwards-compatible, the first type (PEV, BSV or TBV) described
     * (with an "effect" element) will be reported. */
    namedOutMeasures["nMassVaccinations"] =
        OutMeasure::humanDeploy( 22, MHD_VACCINATIONS, Deploy::TIMED );
    /// recoveries in hospital
    namedOutMeasures["nHospitalRecovs"] =
        OutMeasure::humanAC( 23, MHO_HOSPITAL_RECOVERIES, false );
    /// sequelae in hospital
    namedOutMeasures["nHospitalSeqs"] =
        OutMeasure::humanAC( 24, MHO_HOSPITAL_SEQUELAE, false );
    /// Number of IPT Doses (removed together with IPT model)
    namedOutMeasures["nIPTDoses"] = OutMeasure::obsolete( 25 );
    /** Annual Average Kappa
     *
     * Calculated once a year as sum of human infectiousness divided by initial
     * EIR summed over a year. Single value, not per age-group. */
    namedOutMeasures["annAvgK"] = OutMeasure::value( 26, MVF_ANN_AVG_K, true );
    /// Number of episodes (non-malaria fever)
    namedOutMeasures["nNMFever"] =
        OutMeasure::humanAC( 27, MHE_NON_MALARIA_FEVERS, false );
    /// Inoculations per human (all ages) per day of year, over the last year.
    /// (Reporting removed.)
    namedOutMeasures["innoculationsPerDayOfYear"] = OutMeasure::obsolete( 28 );
    /// Kappa (human infectiousness) weighted by availability per day-of-year for the last year.
    /// (Reporting removed.)
    namedOutMeasures["kappaPerDayOfYear"] = OutMeasure::obsolete( 29 );
    /** The total number of inoculations, by age group and cohort, summed over
     * the reporting period. */
    namedOutMeasures["innoculationsPerAgeGroup"] =
        OutMeasure::humanACG( 30, MVF_INOCS, true );
    /// N_v0: emergence of feeding vectors during the last time step. Units: mosquitoes/day
    namedOutMeasures["Vector_Nv0"] = OutMeasure::species( 31, MVF_LAST_NV0, false );
    /// N_v: vectors seeking to feed during the last time step. Units: mosquitoes/day
    namedOutMeasures["Vector_Nv"] = OutMeasure::species( 32, MVF_LAST_NV, false );
    /// N_v: infected vectors seeking to feed during the last time step. Units: mosquitoes/day
    namedOutMeasures["Vector_Ov"] = OutMeasure::species( 33, MVF_LAST_OV, true );
    /// N_v: infectious vectors seeking to feed during the last time step. Units: mosquitoes/day
    namedOutMeasures["Vector_Sv"] = OutMeasure::species( 34, MVF_LAST_SV, true );
    /** Input EIR (Expected EIR entered into scenario file)
     *
     * Units: inoculations per adult per time step. */
    namedOutMeasures["inputEIR"] = OutMeasure::value( 35, MVF_INPUT_EIR, true );
    /** Simulated EIR (EIR output by the transmission model)
     *
     * Units: inoculations per adult per time step (children are excluded
     * when measuring). */
    namedOutMeasures["simulatedEIR"] = OutMeasure::value( 36, MVF_SIM_EIR, true );
    /// Number of Rapid Diagnostic Tests used
    namedOutMeasures["Clinical_RDTs"] = OutMeasure::obsolete( 39 );
    /* Effective total quanty of each drug used orally, in mg.
     * (Per active ingredient abbreviation.)
     * 
     * The quantity is efffective with respect to the cost (see treatment
     * schedule definition).
     * 
     * Reporting removed. */
    namedOutMeasures["Clinical_DrugUsage"] = OutMeasure::obsolete( 40 );
    /// Direct death on first day of CM (before treatment takes effect)
    namedOutMeasures["Clinical_FirstDayDeaths"] =
        OutMeasure::humanAC( 41, MHO_FIRST_DAY_DEATHS, false );
    /// Direct death on first day of CM (before treatment takes effect); hospital only
    namedOutMeasures["Clinical_HospitalFirstDayDeaths"] =
        OutMeasure::humanAC( 42, MHO_HOSPITAL_FIRST_DAY_DEATHS, false );
    /** The number of actual infections since the last survey. */
    namedOutMeasures["nNewInfections"] =
        OutMeasure::humanAC( 43, MHR_NEW_INFECTIONS, false );
    /** The number of ITNs delivered by mass distribution since last survey.
     *
     * These are "modelled ITNs": cover only a single person, cannot be passed
     * to someone else for reuse or used for fishing, etc. */
    namedOutMeasures["nMassITNs"] =
        OutMeasure::humanDeploy( 44, MHD_ITN, Deploy::TIMED );
    /** The number of ITNs delivered through EPI since last survey.
     *
     * Comments from nMassITNs apply. */
    namedOutMeasures["nEPI_ITNs"] =
        OutMeasure::humanDeploy( 45, MHD_ITN, Deploy::CTS );
    /** The number of people newly protected by IRS since last survey.
     *
     * Modelled IRS: affects one person, cannot be plastered over. */
    namedOutMeasures["nMassIRS"] =
        OutMeasure::humanDeploy( 46, MHD_IRS, Deploy::TIMED );
    /** Defunct; was used by "vector availability" intervention (which is now a
     * sub-set of GVI). */
    namedOutMeasures["nMassVA"] = OutMeasure::obsolete( 47 );
    /// Number of malarial tests via microscopy used
    namedOutMeasures["Clinical_Microscopy"] = OutMeasure::obsolete( 48 );
    /* As Clinical_DrugUsage, but for quatities of drug delivered via IV. */
    namedOutMeasures["Clinical_DrugUsageIV"] = OutMeasure::obsolete( 49 );
    /// Number of cohort recruitments removed)
    namedOutMeasures["nAddedToCohort"] = OutMeasure::obsolete( 50 );
    /// Number of individuals removed from cohort (removed)
    namedOutMeasures["nRemovedFromCohort"] = OutMeasure::obsolete( 51 );
    /** Number of people (per age group) treated by mass drug administration
     * campaign. (Note that in one day time-step model MDA can be configured
     * as screen-and-treat. This option reports treatments administered not
     * the number of tests used.) */
    namedOutMeasures["nMDAs"] =
        OutMeasure::humanDeploy( 52, MHD_TREAT, Deploy::TIMED );
    /// Number of deaths caused by non-malaria fevers
    namedOutMeasures["nNmfDeaths"] = OutMeasure::humanAC( 53, MHO_NMF_DEATHS, false );
    /// Number of antibiotic treatments given (disabled — not used)
    namedOutMeasures["nAntibioticTreatments"] = OutMeasure::obsolete( 54 );
    /** Report the number of screenings used in a mass screen-and-treat
     * operation. */
    namedOutMeasures["nMassScreenings"] =
        OutMeasure::humanDeploy( 55, MHD_SCREEN, Deploy::TIMED );
    /// Report the number of mass deployments of generic vector interventions.
    namedOutMeasures["nMassGVI"] =
        OutMeasure::humanDeploy( 56, MHD_GVI, Deploy::TIMED );
    /** Number of IRS deployments via continuous deployment. */
    namedOutMeasures["nCtsIRS"] =
        OutMeasure::humanDeploy( 57, MHD_IRS, Deploy::CTS );
    /** Number of GVI deployments via continuous deployment. */
    namedOutMeasures["nCtsGVI"] =
        OutMeasure::humanDeploy( 58, MHD_GVI, Deploy::CTS );
    /** Number of "MDA" deployments via continuous deployment.
     * 
     * Note: MDA stands for mass drug administration, but the term has come to
     * be used more flexibly by OpenMalaria, including optional screening and
     * deployment through age-based systems. */
    namedOutMeasures["nCtsMDA"] =
        OutMeasure::humanDeploy( 59, MHD_TREAT, Deploy::CTS );
    /** Number of diagnostics used by "MDA" distribution through continuous
     * methods. Can be higher than nCtsMDA since drugs are administered only
     * when the diagnostic is positive. Also see nCtsMDA description. */
    namedOutMeasures["nCtsScreenings"] =
        OutMeasure::humanDeploy( 60, MHD_SCREEN, Deploy::CTS );
    /** Number of removals from a sub-population due to expiry of duration of
     * membership (e.g. intervention too old). */
    namedOutMeasures["nSubPopRemovalTooOld"] =
        OutMeasure::humanAC( 61, MHR_SUB_POP_REM_TOO_OLD, false );
    /** Number of removals from a sub-population due to first
     * infection/bout/treatment (see onFirstBout & co). */
    namedOutMeasures["nSubPopRemovalFirstEvent"] =
        OutMeasure::humanAC( 62, MHR_SUB_POP_REM_FIRST_EVENT, false );
    /** Report the number of Primaquine treatments given. */
    namedOutMeasures["nPQTreatments"] =
        OutMeasure::humanAC( 63, MHT_PQ_TREATMENTS, false );
    /** Report the number of diagnostics used during treatment.
     * 
     * This is not the same as Clinical_RDTs + Clinical_Microscopy: those
     * outputs are used by the "event scheduler" 1-day time step clinical
     * model, whereas this output is used by the 5-day time step model. */
    namedOutMeasures["nTreatDiagnostics"] =
        OutMeasure::humanAC( 64, MHT_TREAT_DIAGNOSTICS, false );
    /** Number of "recruitment only" recruitments via timed deployment. */
    namedOutMeasures["nMassRecruitOnly"] =
        OutMeasure::humanDeploy( 65, MHD_RECRUIT, Deploy::TIMED );
    /** Number of "recruitment only" recruitments via age-based deployment. */
    namedOutMeasures["nCtsRecruitOnly"] =
        OutMeasure::humanDeploy( 66, MHD_RECRUIT, Deploy::CTS );
    /** Number of deployments (of all intervention components) triggered by
     * treatment (case management). */
    namedOutMeasures["nTreatDeployments"] =
        OutMeasure::humanDeploy( 67, MHD_ALL_DEPLOYS, Deploy::TREAT );
    /** Report the total age of all humans in this a group (sum across humans,
     * in years). Divide by nHost to get the average age. */
    namedOutMeasures["sumAge"] = OutMeasure::humanAC( 68, MHF_AGE, true );
    /** The number of human hosts with an infection (patent or not), for each
     * genotype, at the time the survey is taken. */
    namedOutMeasures["nInfectByGenotype"] =
        OutMeasure::humanAC( 69, MHR_INFECTED_GENOTYPE, false );
    /** The number of human hosts whose total (blood-stage) parasite density,
     * for each genotype, is above the detection threshold */
    namedOutMeasures["nPatentByGenotype"] =
        OutMeasure::humanAC( 70, MHR_PATENT_GENOTYPE, false );
}

}
}
