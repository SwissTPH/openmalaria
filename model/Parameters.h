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

#ifndef Hmod_InputParameters
#define Hmod_InputParameters

#include <map>

namespace scnXml { class Parameters; }
namespace OM {

class Parameters {
public:
    enum Parameter {
        /// @b Infection incidence model parameters
        //@{
        NEG_LOG_ONE_MINUS_SINF = 1,
        E_STAR = 2,
        SIMM = 3,
        X_STAR_P = 4,
        GAMMA_P = 5,
        //@}
        /// @b Immunity parameters, mostly on infections
        //@{
        SIGMA_I_SQ = 6,                 ///< Host (not infection) parameter
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
        DENSITY_BIAS_NON_GARKI = 15,        ///< Used in WithinHostModel
        BASELINE_AVAILABILITY_SHAPE = 16,   ///< Used in InfectionIncidenceModel
        LOG_ODDS_RATIO_CF_COMMUNITY = 17,   ///< Used in CaseManagementModel
        INDIRECT_RISK_COFACTOR = 18,        ///< Used in PathogenesisModel
        NON_MALARIA_INFANT_MORTALITY = 19,  ///< Used in Summary
        DENSITY_BIAS_GARKI = 20,            ///< Used in WithinHostModel
        SEVERE_MALARIA_THRESHHOLD = 21,     ///< Used in PathogenesisModel
        IMMUNITY_PENALTY = 22,              ///< Used in WithinHostModel
        IMMUNE_EFFECTOR_DECAY = 23,         ///< Used in WithinHostModel
        /// @b Used in PathogenesisModel
        //@{
        COMORBIDITY_INTERCEPT = 24,
        Y_STAR_HALF_LIFE = 25,
        Y_STAR_1 = 26,
        //@}
        ASEXUAL_IMMUNITY_DECAY = 27,        ///< Used in WithinHostModel
        /// @b Used in PathogenesisModel
        //@{
        Y_STAR_0 = 28,
        
        CRITICAL_AGE_FOR_COMORBIDITY = 30,
        MUELLER_RATE_MULTIPLIER = 31,
        MUELLER_DENSITY_EXPONENT = 32,
        //@}
        /// v in "Case Fatality Rate proposal" TODO: reference
        CFR_SCALE_FACTOR = 33,
        
        // Parameters fitting for Molineaux within host model
        MEAN_LOCAL_MAX_DENSITY = 34,
        SD_LOCAL_MAX_DENSITY = 35,
        MEAN_DIFF_POS_DAYS = 36,
        SD_DIFF_POS_DAYS = 37,
        
        /// exp(-CFR_NEG_LOG_ALPHA) is the proportion of deaths occuring on the first day, with Event Scheduler model
        CFR_NEG_LOG_ALPHA = 38,
        MAX
    };

    Parameters( const scnXml::Parameters& parameters );
    
    /**
     * Get a parameter, using one of the Parameter codes.
     */
    double operator[]( Parameter parameter ) const;
    
private:
    // Initialized (derived) values:
    std::map<Parameter, double> parameterValues;
};

}
#endif
