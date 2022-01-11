/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
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
#include <schema/scenario.h>
#include "util/errors.h"

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
        //@}
        /// @b DescriptiveInfection specific
        //@{
        SIGMA0_SQ = 11,
        X_NU_STAR = 12,
        //@}
        /// @b Used in PathogenesisModel
        //@{
        Y_STAR_SQ = 13,
        ALPHA = 14,
        //@}
        DENSITY_BIAS_NON_GARKI = 15,        ///< Used in Diagnostic
        BASELINE_AVAILABILITY_SHAPE = 16,   ///< Used in InfectionIncidenceModel
        LOG_ODDS_RATIO_CF_COMMUNITY = 17,   ///< Used in CaseManagementModel
        INDIRECT_RISK_COFACTOR = 18,        ///< Used in PathogenesisModel
        NON_MALARIA_INFANT_MORTALITY = 19,  ///< Used in Summary
        DENSITY_BIAS_GARKI = 20,            ///< Used in Diagnostic
        SEVERE_MALARIA_THRESHHOLD = 21,     ///< Used in PathogenesisModel
        IMMUNITY_PENALTY = 22,              ///< Used in WHFalciparum
        IMMUNE_EFFECTOR_DECAY = 23,         ///< Used in WHFalciparum
        /// @b Used in PathogenesisModel
        //@{
        COMORBIDITY_INTERCEPT = 24,
        Y_STAR_HALF_LIFE = 25,
        Y_STAR_1 = 26,
        //@}
        ASEXUAL_IMMUNITY_DECAY = 27,        ///< Used in WHFalciparum
        /// @b Used in PathogenesisModel
        //@{
        Y_STAR_0 = 28,
        
        CRITICAL_AGE_FOR_COMORBIDITY = 30,
        MUELLER_RATE_MULTIPLIER = 31,
        MUELLER_DENSITY_EXPONENT = 32,
        //@}
        /// EventScheduler: v in "Case Fatality Rate proposal"
        CFR_SCALE_FACTOR = 33,
        
        /// @b Molineaux: sampling parameters (not pairwise mode only)
        MEAN_LOCAL_MAX_DENSITY = 34,
        SD_LOCAL_MAX_DENSITY = 35,
        MEAN_DIFF_POS_DAYS = 36,
        SD_DIFF_POS_DAYS = 37,
        
        /// EventScheduler: exp(-CFR_NEG_LOG_ALPHA) is the proportion of deaths occuring on the first day
        CFR_NEG_LOG_ALPHA = 38,
        MAX
    };

    Parameters( const scnXml::Parameters& parameters )
    {
        // set parameters
        const scnXml::Parameters::ParameterSequence& paramSeq = parameters.getParameter();
        for(auto iter = paramSeq.begin(); iter != paramSeq.end(); ++iter) {
            int i = iter->getNumber();
            if (i < 0 || i >= MAX)
                continue;   // ignore the parameter; no real point in making this an error
            Parameter parameter = static_cast<Parameter>(i);
            if( !parameterValues.insert( make_pair( parameter, iter->getValue() ) ).second )
                throw util::xml_scenario_error("parameter with index " + to_string(parameter) + " described twice");
        }
    }

    /**
     * Get a parameter, using one of the Parameter codes.
     */
    double operator[]( Parameter parameter )const{
        auto iter = parameterValues.find( parameter );
        if( iter == parameterValues.end() )
            throw util::xml_scenario_error("parameter " + to_string(parameter) + " required but not described");
        return iter->second;
    }
    
private:
    // Initialized (derived) values:
    std::map<Parameter, double> parameterValues;
};

}
#endif
