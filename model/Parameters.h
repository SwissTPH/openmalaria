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
    /*
    * This defines all possible parameters.
    * Order here is not important.
    */
    enum class Parameter {
        NEG_LOG_ONE_MINUS_SINF,
        E_STAR,
        SIMM,
        X_STAR_P,
        GAMMA_P,
        SIGMA_I_SQ,
        CUMULATIVE_Y_STAR,
        CUMULATIVE_H_STAR,
        NEG_LOG_ONE_MINUS_ALPHA_M,
        DECAY_M,
        SIGMA0_SQ,
        X_NU_STAR,
        Y_STAR_SQ,
        ALPHA,
        DENSITY_BIAS_NON_GARKI,
        BASELINE_AVAILABILITY_SHAPE,
        LOG_ODDS_RATIO_CF_COMMUNITY,
        INDIRECT_RISK_COFACTOR,
        NON_MALARIA_INFANT_MORTALITY,
        DENSITY_BIAS_GARKI,
        SEVERE_MALARIA_THRESHHOLD,
        IMMUNITY_PENALTY,
        IMMUNE_EFFECTOR_DECAY,
        COMORBIDITY_INTERCEPT,
        Y_STAR_HALF_LIFE,
        Y_STAR_1,
        ASEXUAL_IMMUNITY_DECAY,
        Y_STAR_0,
        CRITICAL_AGE_FOR_COMORBIDITY,
        MUELLER_RATE_MULTIPLIER,
        MUELLER_DENSITY_EXPONENT,
        CFR_SCALE_FACTOR,
        MEAN_LOCAL_MAX_DENSITY,
        SD_LOCAL_MAX_DENSITY,
        MEAN_DIFF_POS_DAYS,
        SD_DIFF_POS_DAYS,
        CFR_NEG_LOG_ALPHA
    };

    /*
    * Initializes parameters using what is specified in the input XML, and does some
    * validation on the specified parameters and values.
    */
    Parameters( const scnXml::Parameters& parameters )
    {
        const scnXml::Parameters::ParameterSequence& paramSeq = parameters.getParameter();
        for(auto iter = paramSeq.begin(); iter != paramSeq.end(); ++iter) {
            const int paramId = iter->getNumber();
            const double paramValue = iter->getValue();

            // C++20 offers "contains" which would be more expressive than "count" here.
            const bool paramIdIsValid = IdToParamMap.count(paramID);
            const bool paramIsAlreadyDefined = parameterValues.count(paramID);
            if (!paramIdIsValid)
            {
                throw util::xml_scenario_error("index " + to_string(paramID) +
                        " does not correspond to any parameter")
            }
            else if (paramIsAlreadyDefined)
            {
                throw util::xml_scenario_error("parameter with index " + to_string(parameter) +
                        " described twice");
            }
            else
            {
                Parameter targetParam = IdToParamMap[paramID];
                paramValues[targetParam] = parameterValue;
            }
        }
    }

    /**
     * Get a parameter, using one of the Parameter codes.
     */
    double operator[]( Parameter parameter )const{
        auto iter = parameterValues.find( parameter );
        if( iter == parameterValues.end() )
            throw util::xml_scenario_error("parameter " + to_string(parameter) +
                    " required but not described");
        return iter->second;
    }

private:
    // Initialized (derived) values:
    std::map<Parameter, double> parameterValues;

    /*
    * When a user sets a parameter's value in XML, the important things for them to specify are
    * the integer ID of the parameter, and the value they want the parameter to have.
    * This map defines the relation between such integer IDs and actual parameters.
    */
    std::map<int, Parameter> IdToParamMap {
        /// @b Infection incidence model parameters
        //@{
        { 1, NEG_LOG_ONE_MINUS_SINF },
        { 2, E_STAR },
        { 3, SIMM },
        { 4, X_STAR_P },
        { 5, GAMMA_P },

        //@}
        /// @b Immunity parameters, mostly on infections
        //@{
        { 6, SIGMA_I_SQ },                 ///< Host (not infection) parameter
        { 7, CUMULATIVE_Y_STAR },
        { 8, CUMULATIVE_H_STAR },
        { 9, NEG_LOG_ONE_MINUS_ALPHA_M },
        { 10, DECAY_M },

        //@}
        /// @b DescriptiveInfection specific
        //@{
        { 11, SIGMA0_SQ },
        { 12, X_NU_STAR },

        //@}
        /// @b Used in PathogenesisModel
        //@{
        { 13, Y_STAR_SQ },
        { 14, ALPHA },
        //@}

        { 15, DENSITY_BIAS_NON_GARKI },        ///< Used in Diagnostic
        { 16, BASELINE_AVAILABILITY_SHAPE },   ///< Used in InfectionIncidenceModel
        { 17, LOG_ODDS_RATIO_CF_COMMUNITY },   ///< Used in CaseManagementModel
        { 18, INDIRECT_RISK_COFACTOR },        ///< Used in PathogenesisModel
        { 19, NON_MALARIA_INFANT_MORTALITY },  ///< Used in Summary
        { 20, DENSITY_BIAS_GARKI },            ///< Used in Diagnostic
        { 21, SEVERE_MALARIA_THRESHHOLD },     ///< Used in PathogenesisModel
        { 22, IMMUNITY_PENALTY },              ///< Used in WHFalciparum
        { 23, IMMUNE_EFFECTOR_DECAY },         ///< Used in WHFalciparum

        /// @b Used in PathogenesisModel
        //@{
        { 24, COMORBIDITY_INTERCEPT },
        { 25, Y_STAR_HALF_LIFE },
        { 26, Y_STAR_1 },
        //@}

        { 27, ASEXUAL_IMMUNITY_DECAY },        ///< Used in WHFalciparum

        /// @b Used in PathogenesisModel
        //@{
        { 28, Y_STAR_0 },
        { 30, CRITICAL_AGE_FOR_COMORBIDITY },
        { 31, MUELLER_RATE_MULTIPLIER },
        { 32, MUELLER_DENSITY_EXPONENT },
        //@}

        /// EventScheduler: v in "Case Fatality Rate proposal"
        { 33, CFR_SCALE_FACTOR },

        /// @b Molineaux: sampling parameters (not pairwise mode only)
        { 34, MEAN_LOCAL_MAX_DENSITY },
        { 35, SD_LOCAL_MAX_DENSITY },
        { 36, MEAN_DIFF_POS_DAYS },
        { 37, SD_DIFF_POS_DAYS },

        /// EventScheduler: exp(-CFR_NEG_LOG_ALPHA) is the proportion of deaths occuring on the first day
        { 38, CFR_NEG_LOG_ALPHA }
    };
};

}
#endif
