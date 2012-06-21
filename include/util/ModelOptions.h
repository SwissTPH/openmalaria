/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef Hmod_util_ModelOptions
#define Hmod_util_ModelOptions

#include "Global.h"

class UnittestUtil;

namespace OM { namespace util {
    
    /** Flags signalling which versions of some models to use. */
    enum OptionCodes {
	/** @brief Clinical episodes reduce the level of acquired immunity
	* 
	* Effective cumulative exposure to blood stage parasites is reduced by a
	* clinical sickness event, so that clinical bouts have a negative effect on
	* blood stage immunity.
	* 
	* (ImmediateOutcomes model: per event; EventScheduler: once per bout.)
	* 
	* Default: Clinical events have no effect on immune status except
	* secondarily via effects of treatment. */
	PENALISATION_EPISODES = 0,
	
	/** @brief Baseline availability of humans is sampled from a gamma distribution
	* Infections introduced by mass action with negative binomial
	* variation in numbers of infection.
	* 
	* Default: New infections are introduced via a Poisson process as described
	* in AJTMH 75 (suppl 2) pp11-18. */
	NEGATIVE_BINOMIAL_MASS_ACTION,
	
	/** @brief An IPT model, no longer used
	* 
	* Does nothing if IPT is not present. */
	ATTENUATION_ASEXUAL_DENSITY,
	
	/** @brief Baseline availability of humans is sampled from a log normal distribution
	* 
	* Infections introduced by mass action with log normal variation in
	* infection rate.
	* 
	* Default: New infections are introduced via a Poisson process as described
	* in AJTMH 75 (suppl 2) pp11-18. */
	LOGNORMAL_MASS_ACTION,
	
	/** Infections are introduced without using preerythrocytic immunity. */
	NO_PRE_ERYTHROCYTIC,
	
	/** @brief Bug fixes in Descriptive & DescriptiveIPT within-host models.
	*
	* Really, the code makes more sense if all these are used. But, they are
	* not always used in order to preserve consistant results.
	* 
	* MAX_DENS_RESET is not used since it is unneeded when MAX_DENS_CORRECTION
	* is present and wouldn't make sense when not. */
	// @{
	MAX_DENS_CORRECTION,
	INNATE_MAX_DENS,
// 	MAX_DENS_RESET,
	//@}
	
	/** @brief Parasite densities are predicted from an autoregressive process
	*
	* Default: Parasite densities are determined from the descriptive model
	* given in AJTMH 75 (suppl 2) pp19-31 .*/
	DUMMY_WITHIN_HOST_MODEL,
	
	/** Clinical episodes occur if parasitaemia exceeds the pyrogenic threshold.
	* 
	* Default: Clinical episodes are a stochastic function as described in AJTMH
	* 75 (suppl 2) pp56-62. */
	PREDETERMINED_EPISODES,
	
	/** @brief The presentation model includes simulation of non-malaria fevers
	* 
	* Default: Non-malaria fevers are not simulated. */
	NON_MALARIA_FEVERS,
	
	/** @brief Use a PK & PD model for drug effects
	 *
	 * This causes simulation of the pharmacokinetics and pharmacodynamics of drugs,
	 * as opposed to the original (and default) models, in which drugs have all or
	 * nothing effects (except in certain IPTi models).
	 * 
	 * Currently this means to use either the Hoshen or the LSTM PKPD model
	 * (LSTM when drugDescription XML element is present in scenario description.)
	 */
	INCLUDES_PK_PD,
	
	/** @brief Use revised clinical and case management model, ClinicalEventScheduler
	* 
	* Default: use the Tediosi et al case management model (Case management as
	* described in AJTMH 75 (suppl 2) pp90-103), ClinicalImmediateOutcomes. */
	CLINICAL_EVENT_SCHEDULER,
	
	/** @brief Clinical episodes occur in response to a simple parasite density trigger
	* 
	* Default: Use the Ross et al presentation model (Clinical episodes are a
	* stochastic function as described in AJTMH 75 (suppl 2) pp56-62). */
	MUELLER_PRESENTATION_MODEL,
	
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
	TRANS_HET,
	/// @brief Allow simple heterogeneity in comorbidity
	COMORB_HET,
	/// @brief Allow simple heterogeneity in treatment seeking
	TREAT_HET,
	/// @brief Allow correlated heterogeneities in transmission and comorbidity
	COMORB_TRANS_HET,
	/// @brief Allow correlated heterogeneities in transmission and treatment seeking
	TRANS_TREAT_HET,
	/// @brief Allow correlated heterogeneities comorbidity and treatment seeking
	COMORB_TREAT_HET,
	/// @brief Allow correlated heterogeneities in transmission, comorbidity and treatment seeking
	TRIPLE_HET,
        // @}
        
	/** @brief Selection of within host models
	*/
        //@{
        /// Parasite densities are predicted from an empirical model
	EMPIRICAL_WITHIN_HOST_MODEL,
	
        /// Use Molineaux within host model
	MOLINEAUX_WITHIN_HOST_MODEL,
        
        /** Use Penny infection model. */
        PENNY_WITHIN_HOST_MODEL,
        // @}
	
	/** @brief Selection of gamma distribution for Molineaux or Penny within host models
	 */
	//@{
	/** Use Gamma distribution for mean duration (Molineaux model) */
	MEAN_DURATION_GAMMA,
	/** Use Gamma distribution for first local maximum (Molineaux model) */
	FIRST_LOCAL_MAXIMUM_GAMMA,
	/** Use Gamma distribution for first local maximum (Molineaux model) */
	PARASITE_REPLICATION_GAMMA,
	
	/** Use Gamma distribution for immune threshold (Penny model) */
	IMMUNE_THRESHOLD_GAMMA,
	
	/** Use Gamma distribution for update density (Penny model) */
	UPDATE_DENSITY_GAMMA,
	//@}
        
	/** Use the Garki density bias instead of the default one in the detection limit.
	 *
	 * The default bias corresponds to counting parasites and white blood cells
	 * (assuming a white blood cell density of 8000 per µl), the Garki bias to
	 * estimations from a probability function. */
	GARKI_DENSITY_BIAS,
	
	/** Use the IPT(i) drug model (DescriptiveIPTWithinHost and
	 * DescriptiveIPTInfection classes) with its simple SP model.
	 * 
	 * NOTE: code is unmaintained in order to keep results comparable with
	 * previous experiments run. On a 1-day timestep, we now have a better
	 * drug model, but will need to rewrite IPT code to purely be an
	 * intervention.
         * 
         * Note: previously this implied REPORT_ONLY_AT_RISK; now it does not.
	 */
	IPTI_SP_MODEL,
	
	/** Turn off reporting of several outputs for humans suffering a recent
         * clinical episode and therefore not currently at risk of what would
         * clinically be regarded as a separate episode.
         * 
         * <b>This is a compatibility option only.</b> It only works with the
         * 5-day model and the length of the not-at-risk period is hard-coded,
         * not dependant on the healthSystemMemory value.
         * 
         * This removes from several outputs all humans who recieved treatment
         * during the previous 4 (5-day) timesteps, who are therefore not
         * currently at risk of an additional episode. Summaries affected include
         * nHost, nInfect, nExpectd, nPatent, totalInfs, totalPatentInf, sumlogDens,
         * nNewInfections, sumLogPyrogenThres, sumPyrogenThresh, and potentially
         * other outputs added after writing this. */
	REPORT_ONLY_AT_RISK,
        
        /** Turn on vector life-cycle model. Allows better larviciding and
         * vector population dynamics modelling.
         * 
         * Requires vector model. */
        VECTOR_LIFE_CYCLE_MODEL,
        
        /** Turn on the simple mosquito population dynamics model (simpler
         * version of life-cycle model).
         * 
         * Requires vector model. */
        VECTOR_SIMPLE_MPD_MODEL,
        
	// Used by tests; should be 1 more than largest option
	NUM_OPTIONS,
    };
    
    
    /// Encapsulation for "modelVersion" xml attribute
    class ModelOptions {
    public:
	/** Return true if given option (from OptionCodes) is active. */
	static inline bool option(OptionCodes code) {
	    /* Performance:
	    using bitset<>::operator[] constructs and returns some kind of
	    reference. It's slow! bitset<>::test() isn't much better, hence
	    reverted to integer binary ops (1/8th time). (Also note: this was
	    the only use of bitset<> with significant performance impact.)
	    */
	    return (optSet & (1<<code)) != 0u;
	}
	/** Return true if any of TRANS_HET, COMORB_TRANS_HET, TRANS_TREAT_HET or
	 * TRIPLE_HET are active. */
	static inline bool anyTransHet () {
	    static const uint32_t anyHet = 1<<TRANS_HET | 1<<COMORB_TRANS_HET | 1<<TRANS_TREAT_HET | 1<<TRIPLE_HET;
	    return (optSet & anyHet) != 0u;
	}
	
	/// Set options from XML file.
        /// Relies on Global::init() already having been called.
	static void init ();
	
    private:
	/** Model options.
	 *
	 * Default value set by init(). */
	static uint32_t optSet;
	
	friend class ::UnittestUtil;	// Note: class is in base namespace
    };
} }
#endif
