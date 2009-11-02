/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef CONSTANT_H
#define CONSTANT_H

/** Flags signalling which versions of some models to use. */
enum ModelVersion {
  /* Values are written here using left-shifts. 1 << x is equal to pow(2,x)
   * for integers, so each value here has only one bit true in binary, allowing
   * the bits to be used as flags: http://en.wikipedia.org/wiki/Flag_byte
   * For historical reasons only, there is no 1 (=1<<0). */
  /** @brief Clinical episodes reduce the level of acquired immunity
   * 
   * Effective cumulative exposure to blood stage parasites is reduced during a
   * clinical episode, so that clinical episodes have a negative effect on
   * blood stage immunity.
   * 
   * Default: Clinical events have no effect on immune status except
   * secondarily via effects of treatment. */
  PENALISATION_EPISODES = 1 << 1,
  
  /** @brief Baseline availability of humans is sampled from a gamma distribution
   * Infections introduced by mass action with negative binomial
   * variation in numbers of infection.
   * 
   * Default: New infections are introduced via a Poisson process as described
   * in AJTMH 75 (suppl 2) pp11-18. */
  NEGATIVE_BINOMIAL_MASS_ACTION = 1 << 2,
  
  /** @brief 
   * 
   * Does nothing if IPT is not present. */
  ATTENUATION_ASEXUAL_DENSITY = 1 << 3,
  
  /** @brief Baseline availability of humans is sampled from a log normal distribution
   * 
   * Infections introduced by mass action with log normal variation in
   * infection rate.
   * 
   * Default: New infections are introduced via a Poisson process as described
   * in AJTMH 75 (suppl 2) pp11-18. */
  LOGNORMAL_MASS_ACTION = 1 << 4,
  
  /** Infections are introduced without using preerythrocytic immunity. */
  NO_PRE_ERYTHROCYTIC = 1 << 5,
  
  /// BugFix in previous versions.  This option is not currently implemented.
  // @{
  MAX_DENS_CORRECTION = 1 << 6,
  INNATE_MAX_DENS = 1 << 7,
  MAX_DENS_RESET = 1 << 8,
  //@}
  
  /** @brief Parasite densities are predicted from an autoregressive process
   *
   * Default: Parasite densities are determined from the descriptive model
   * given in AJTMH 75 (suppl 2) pp19-31 .*/
  DUMMY_WITHIN_HOST_MODEL = 1 << 9,
  
  /** Clinical episodes occur if parasitaemia exceeds the pyrogenic threshold.
   * 
   * Default: Clinical episodes are a stochastic function as described in AJTMH
   * 75 (suppl 2) pp56-62. */
  PREDETERMINED_EPISODES = 1 << 10,
  
  /** @brief The presentation model includes simulation of non-malaria fevers
   * 
   * Default: Non-malaria fevers are not simulated. */
  NON_MALARIA_FEVERS = 1 << 11,
  
  /** @brief Pharmacokinetic and pharmacodynamics of drugs are simulated
   * 
   * Default: Drugs have all or nothing effects (except in certain IPTi
   * models). */
  INCLUDES_PK_PD = 1 << 12,
  
  /** @brief Use revised clinical and case management model, ClinicalEventScheduler
   * 
   * Default: use the Tediosi et al case management model (Case management as
   * described in AJTMH 75 (suppl 2) pp90-103), ClinicalImmediateOutcomes. */
  CLINICAL_EVENT_SCHEDULER = 1 << 13,
  
  /** @brief Clinical episodes occur in response to a simple parasite density trigger
   * 
   * Default: Use the Ross et al presentation model (Clinical episodes are a
   * stochastic function as described in AJTMH 75 (suppl 2) pp56-62). */
  MUELLER_PRESENTATION_MODEL = 1 << 14,
  
  /** @brief Simple heterogeneity
   * 
   * Defaults: No heterogeneity.
   * 
   * (Transmission) heterogeneity is incompatible with
   * NEGATIVE_BINOMIAL_MASS_ACTION and LOGNORMAL_MASS_ACTION because both try
   * to adjust _EIRFactor and it is not confirmed that the ways they do this is
   * compatible. */
  // @{
  /// @brief Allow simple heterogeneity in transmission
  TRANS_HET = 1 << 15,
  /// @brief Allow simple heterogeneity in comorbidity
  COMORB_HET = 1 << 16,
  /// @brief Allow simple heterogeneity in treatment seeking
  TREAT_HET = 1 << 17,
  /// @brief Allow correlated heterogeneities in transmission and comorbidity
  COMORB_TRANS_HET = 1 << 18,
  /// @brief Allow correlated heterogeneities in transmission and treatment seeking
  TRANS_TREAT_HET = 1 << 19,
  /// @brief Allow correlated heterogeneities comorbidity and treatment seeking
  COMORB_TREAT_HET = 1 << 20,
  /// @brief Allow correlated heterogeneities in transmission, comorbidity and treatment seeking
  TRIPLE_HET = 1 << 21,

  /** @brief Parasite densities are predicted from an empirical model
   */
  EMPIRICAL_WITHIN_HOST_MODEL = 1 << 22,

  /// Used to test if any heterogeneity is present
  ANY_HET = TRANS_HET|COMORB_HET|TREAT_HET|COMORB_TRANS_HET|TRANS_TREAT_HET|TRIPLE_HET,
  ANY_TRANS_HET =  TRANS_HET | COMORB_TRANS_HET | TRANS_TREAT_HET | TRIPLE_HET,
  // @}
  
  // Used by tests; should be 1 plus highest left-shift value of 1
  NUM_VERSIONS = 23,
};

/// Namespace enclosing pathogenesis output enumeration.
namespace Pathogenesis {
  /** Types of sickness; used by case management.
   *
   * Most values are flags which can be combined in any form. A few
   * combinations set follow. */
  enum State {
    /* Values here are written in hexadecimal: http://en.wikipedia.org/wiki/Hexadecimal
     * Many are designed to be "flags", so the value corresponds to a single bit:
     * http://en.wikipedia.org/wiki/Flag_byte
     * Max: 0x4000
     * (note & | ^ are C++'s binary AND, OR and XOR operators). */
    NONE		= 0,		///< Not sick
    
    // Flags for current state/worst state to report:
    SICK		= 0x1,		///< Sick (may or may not be from malaria)
    MALARIA		= 0x2,		///< Malaria sickness
    SEVERE		= 0x8,		///< Severe malaria case
    COINFECTION		= 0x4,		///< Malaria with a coinfection
    /// Used by ClinicalEventScheduler to indicate a second bout of malarial sickness within the same episode (roughly)
    SECOND_CASE		= 0x10,
    COMPLICATED		= 0x20,		///< Flag used to indicate SEVERE and/or COINFECTION
    
    // Flag used by pathogenesis model to tell the clinical model that individual will die; not used for reporting:
    INDIRECT_MORTALITY	= 0x800,	///< Death caused by indirect effects of malaria
    
    // Flags for outcome reporting:
    EVENT_IN_HOSPITAL	= 0x400,	///< Indicates recovery/sequelae/death event occurred in hospital − only set on one of these events
    DIRECT_DEATH	= 0x1000,	///< Used for reporting death (from COMPLICATED sickness)
    SEQUELAE		= 0x2000,	///< Reporting recovered with sequelae (from COMPLICATED sickness)
    RECOVERY		= 0x4000,	///< Report that individual fully recovered
    
    STATE_MALARIA	= SICK | MALARIA,	///< Combination: SICK, MALARIA
    STATE_SEVERE	= STATE_MALARIA | COMPLICATED | SEVERE,	///< Combination: SICK, MALARIA, COMPLICATED, SEVERE
    STATE_COINFECTION	= STATE_MALARIA | COMPLICATED | COINFECTION,	///< Combination: SICK, MALARIA, COMPLICATED, COINFECTION
  };
}

namespace Params {
  enum Params {
    /// @b Used in NoVectorControl
    //@{
    NEG_LOG_ONE_MINUS_SINF = 1,
    E_STAR = 2,
    SIMM = 3,
    X_STAR_P = 4,
    GAMMA_P = 5,
    //@}
    SIGMA_I_SQ = 6,			///< Used in WithinHostModel
    /// @b Used in Infection
    //@{
    CUMULATIVE_Y_STAR = 7,
    CUMULATIVE_H_STAR = 8,
    NEG_LOG_ONE_MINUS_ALPHA_M = 9,
    DECAY_M = 10,
    SIGMA0_SQ = 11,
    X_NU_STAR = 12,
    //@}
    /// @b Used in PathogenesisModel
    //@{
    Y_STAR_SQ = 13,
    ALPHA = 14,
    //@}
    DENSITY_BIAS_NON_GARKI = 15,	///< Used in WithinHostModel
    BASELINE_AVAILABILITY_SHAPE = 16,	///< Used in InfectionIncidenceModel
    LOG_ODDS_RATIO_CF_COMMUNITY = 17,	///< Used in CaseManagementModel
    INDIRECT_RISK_COFACTOR = 18,	///< Used in PathogenesisModel
    NON_MALARIA_INFANT_MORTALITY = 19,	///< Used in Summary
    DENSITY_BIAS_GARKI = 20,		///< Used in WithinHostModel
    SEVERE_MALARIA_THRESHHOLD = 21,	///< Used in PathogenesisModel
    IMMUNITY_PENALTY = 22,		///< Used in WithinHostModel
    IMMUNE_EFFECTOR_DECAY = 23,		///< Used in WithinHostModel
    /// @b Used in PathogenesisModel
    //@{
    COMORBIDITY_INTERCEPT = 24,
    Y_STAR_HALF_LIFE = 25,
    Y_STAR_1 = 26,
    //@}
    ASEXUAL_IMMUNITY_DECAY = 27,	///< Used in WithinHostModel
    /// @b Used in PathogenesisModel
    //@{
    Y_STAR_0 = 28,
    
    CRITICAL_AGE_FOR_COMORBIDITY = 30,
    MUELLER_RATE_MULTIPLIER = 31,
    MUELLER_DENSITY_EXPONENT = 32,
    //@}
    MAX
  };
}

/** Value used as the timestep for an event which has never happened.
 *
 * For any simulation timestep, we must have:
 * ( TIMESTEP_NEVER + simulationTime < 0 )
 * but since (x - TIMESTEP_NEVER >= y) is often checked, x - TIMESTEP_NEVER
 * must not overflow for any timestep x (int represents down to -0x7FFFFFFF).
 */
const int TIMESTEP_NEVER = -0x3FFFFFFF;

/// Days in a year. Should be a multiple of interval.
const int daysInYear= 365;

/** There are 3 simulation modes. */
enum SimulationMode {
  /** Equilibrium mode
   * 
   * This is used for the warm-up period and if we want to separate direct
   * effect of an intervention from indirect effects via transmission
   * intensity. The seasonal pattern and intensity of the EIR do not change
   * over years.
   * 
   * For the vector model, this runs most calculations dynamically but still
   * forces the EIR. */
  equilibriumMode = 2,
  
  /** Transient EIR known
   * 
   * This is used to simulate an intervention that changes EIR, and where we
   * have measurements of the EIR over time during the intervention period. */
  transientEIRknown = 3,
  
  /** EIR changes
   * 
   * The simulation is driven by the EIR which changes dynamically during the
   * intervention phase as a function of the characteristics of the
   * interventions.
   * 
   * Dependending on whether the Vector or NonVector model is in use, this EIR
   * may be calculated from a mosquito emergence rate or be an input EIR
   * scaled by the relative infectiousness of the humans. */
  dynamicEIR = 4,
};

#endif
