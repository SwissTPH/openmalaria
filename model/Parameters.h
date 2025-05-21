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
#include <optional>
#include <unordered_map>
#include <schema/scenario.h>
#include "util/CommandLine.h"
#include "util/errors.h"

namespace scnXml { class Parameters; }
namespace OM {

/*
* Defines the names that clients of the Parameters class use in order to
* read the values of parameters.
*
* NOTE: any time that a new parameter is added here, a corresponding entry
* *must* be added to the map in the Parameters class.  Such an entry relates
* a parameter codes (positive ingeger) to the new parameter names.  If no
* such entry is added, the new parameter will not be usable in the simulation.
*/
enum class ParameterName {
    /// @b Infection incidence model parameters
    //@{
    NEG_LOG_ONE_MINUS_SINF,
    E_STAR,
    SIMM,
    X_STAR_P,
    GAMMA_P,
    //@}
    /// @b Immunity parameters, mostly on infections
    //@{
    SIGMA_I_SQ,                 ///< Host (not infection) parameter
    CUMULATIVE_Y_STAR,
    CUMULATIVE_H_STAR,
    NEG_LOG_ONE_MINUS_ALPHA_M,
    DECAY_M,
    //@}
    /// @b DescriptiveInfection specific
    //@{
    SIGMA0_SQ,
    X_NU_STAR,
    //@}
    /// @b Used in PathogenesisModel
    //@{
    Y_STAR_SQ,
    ALPHA,
    //@}
    DENSITY_BIAS_NON_GARKI,        ///< Used in Diagnostic
    BASELINE_AVAILABILITY_SHAPE,   ///< Used in InfectionIncidenceModel
    LOG_ODDS_RATIO_CF_COMMUNITY,   ///< Used in CaseManagementModel
    INDIRECT_RISK_COFACTOR,        ///< Used in PathogenesisModel
    NON_MALARIA_INFANT_MORTALITY,  ///< Used in Summary
    DENSITY_BIAS_GARKI,            ///< Used in Diagnostic
    SEVERE_MALARIA_THRESHHOLD,     ///< Used in PathogenesisModel
    IMMUNITY_PENALTY,              ///< Used in WHFalciparum
    IMMUNE_EFFECTOR_DECAY,         ///< Used in WHFalciparum
    /// @b Used in PathogenesisModel
    //@{
    COMORBIDITY_INTERCEPT,
    Y_STAR_HALF_LIFE,
    Y_STAR_1,
    //@}
    ASEXUAL_IMMUNITY_DECAY,        ///< Used in WHFalciparum
    /// @b Used in PathogenesisModel
    //@{
    Y_STAR_0,

    CRITICAL_AGE_FOR_COMORBIDITY,
    MUELLER_RATE_MULTIPLIER,
    MUELLER_DENSITY_EXPONENT,
    //@}
    /// EventScheduler: v in "Case Fatality Rate proposal"
    CFR_SCALE_FACTOR,

    /// @b Molineaux: sampling parameters (not pairwise mode only)
    MEAN_LOCAL_MAX_DENSITY,
    SD_LOCAL_MAX_DENSITY,
    MEAN_DIFF_POS_DAYS,
    SD_DIFF_POS_DAYS,

    /// EventScheduler: exp(-CFR_NEG_LOG_ALPHA) is the proportion of deaths occuring on the first day
    CFR_NEG_LOG_ALPHA,
};

/*
* This class' job is to make the values of model parameters available to clients when the simulation
* is running.
*
* From a user perspective, each parameter has these 3 parts:
*
* 1. The name of the parameter as written in XML.
*       OpenMalaria completely disregards this.
*
* 2. The numerical ID, "index" which identfies the parameter in XML.
*       OpenMalaria relies on this.
*
* 3. The value assigned to the parameter in XML.
*       OpenMalaria associates this with the numerical ID.
*
* From a developer perspective, each "parameter" can be thought of as having these 3 parts:
*
* 1. A ParameterName.
*       Clients of this class use a ParameterName to refer to the parameter.
*
* 2. A numerical ID to identify the parameter.
*       This is used to map parameters in XML to parameters in OpenMalaria,
*       and for error reporting.
*
* 2. A (possible) floating point number representing the actual value assigned to the parameter. 
*       If the user does not explicitly or implicitly set a value, an empty
*       value is assigned.  If a client of this class attempts to read the
*       parameter's value, but the value is empty, this is considered an
*       error - and this class is responsible for handling such errors.
*/
class Parameters {
public:

    Parameters( const scnXml::Parameters& parameters )
    {
        for (const std::pair<int, ParameterName> pair : idCodeToNameMap)
        {
            const ParameterName name = pair.second;
            nameToValueMap[name] = std::nullopt;
        }

        initializeParamsFromXML(parameters);
    }

    /**
     * Get a parameter, using one of the Parameter codes.
     */
    double operator[]( ParameterName name )const{
        // First check the parameter actually exists and has a value.
        const bool paramHasAValueSet = nameToValueMap.at(name) != std::nullopt;
        if (!paramHasAValueSet)
        {
            for (const std::pair<int, ParameterName> pair : idCodeToNameMap)
            {
                if (pair.second == name)
                {
                    const int paramCode = pair.first;
                    throw util::xml_scenario_error("Parameter with index " + to_string(paramCode) +
                            " required but not described.");
                }
            }
            throw util::base_exception("A parameter required by this simulation is missing a definition fo r its ID in C++ code.");
        }

        return nameToValueMap.at(name).value();
    }

private:

    /*
    * Initializes parameters using explicit values specified in the input XML,
    * and does some validation on the specified parameters and values.
    * This method will be used if the user does not use the "base" model or any
    * other pre-set collection of parameter values.
    */
    void initializeParamsFromXML(const scnXml::Parameters& parameters)
    {
        const scnXml::Parameters::ParameterSequence& paramSeq = parameters.getParameter();

        for(auto iter = paramSeq.begin(); iter != paramSeq.end(); ++iter)
        {
            const int paramId = iter->getNumber();
            const double paramValueFromXML = iter->getValue();

            // C++20 offers "contains" which would be more expressive than "count" here.
            const bool paramIdIsValid = idCodeToNameMap.count(paramId);
            if (!paramIdIsValid)
            {
                if( util::CommandLine::option(util::CommandLine::DEPRECATION_WARNINGS) ){
                    cerr << "Deprecation warning: <parameter> index" << to_string(paramId) <<
                        " is no longer used" << endl;
                }

                // If we throw here, old scenarios containing the deprecated parameter
                // will need manual work to migrate to new OpenMalaria versions.  This
                // isn't unacceptable, but it's simpler here to just skip over the
                // deprecated parameter.  In particular, many tests specify values for
                // deprecated some parameters.
                continue;
            }
            const ParameterName nameOfParamToSet = idCodeToNameMap.at(paramId);

            // C++20 offers "contains" which would be more expressive than "count" here.
            const bool paramValueAlreadySet = nameToValueMap.at(nameOfParamToSet) != std::nullopt;
            if (paramValueAlreadySet)
            {
                throw util::xml_scenario_error("Parameter with index " + to_string(paramId) +
                        " described twice.");
            }

            nameToValueMap[nameOfParamToSet] = paramValueFromXML;
        }
    }

    /*
    * Initializes some hardcoded values for some parameters, encapsulating the parameter values
    * that make up the base model.
    */
    void initializeParamsBaseModel()
    {
           nameToValueMap[ idCodeToNameMap.at( 1) ] = 0.050736; // '-ln(1-Sinf)'
           nameToValueMap[ idCodeToNameMap.at( 2) ] = 0.03247; // Estar
           nameToValueMap[ idCodeToNameMap.at( 3) ] = 0.138161050830301; // Simm
           nameToValueMap[ idCodeToNameMap.at( 4) ] = 1514.385853233699891; // Xstar_p
           nameToValueMap[ idCodeToNameMap.at( 5) ] = 2.03692533424484; // gamma_p
           nameToValueMap[ idCodeToNameMap.at( 6) ] = 10.173598698525799; // sigma2i
           nameToValueMap[ idCodeToNameMap.at( 7) ] = 35158523.31132510304451; // CumulativeYstar
           nameToValueMap[ idCodeToNameMap.at( 8) ] = 97.334652723897705; // CumulativeHstar
           nameToValueMap[ idCodeToNameMap.at( 9) ] = 2.33031045876193; // '-ln(1-alpha_m)'
           nameToValueMap[ idCodeToNameMap.at(10) ] = 2.53106547375805; // decay_m
           nameToValueMap[ idCodeToNameMap.at(11) ] = 0.655747311168152; // sigma2_0
           nameToValueMap[ idCodeToNameMap.at(12) ] = 0.916181104713054; // Xstar_v
           nameToValueMap[ idCodeToNameMap.at(13) ] = 6502.26335600001039; // Ystar2
           nameToValueMap[ idCodeToNameMap.at(14) ] = 142601.912520000012591; // alpha
           nameToValueMap[ idCodeToNameMap.at(15) ] = 0.177378570987455; // Density bias (non Garki)
           nameToValueMap[ idCodeToNameMap.at(16) ] = 1.0; //  sigma2
           nameToValueMap[ idCodeToNameMap.at(17) ] = 0.736202; // log oddsr CF community
           nameToValueMap[ idCodeToNameMap.at(18) ] = 0.018777338; // Indirect risk cofactor
           nameToValueMap[ idCodeToNameMap.at(19) ] = 49.539046599999999; // Non-malaria infant mortality
           nameToValueMap[ idCodeToNameMap.at(20) ] = 4.79610772546704; // Density bias (Garki)
           nameToValueMap[ idCodeToNameMap.at(21) ] = 784455.599999999976717; // Severe Malaria Threshhold
           nameToValueMap[ idCodeToNameMap.at(22) ] = 1; // Immunity Penalty
           nameToValueMap[ idCodeToNameMap.at(23) ] = 0; // Immune effector decay
           nameToValueMap[ idCodeToNameMap.at(24) ] = 0.0968; // comorbidity intercept
           nameToValueMap[ idCodeToNameMap.at(25) ] = 0.275437402; // Ystar half life
           nameToValueMap[ idCodeToNameMap.at(26) ] = 0.596539864; // Ystar1
           nameToValueMap[ idCodeToNameMap.at(27) ] = 0; // Asexual immunity decay
           nameToValueMap[ idCodeToNameMap.at(28) ] = 296.302437899999973; // Ystar0
           nameToValueMap[ idCodeToNameMap.at(29) ] = 2.797523626; // Idete multiplier
           nameToValueMap[ idCodeToNameMap.at(30) ] = 0.117383; // critical age for comorbidity
    }

    /*
    * Defines the map from parameter ID numbers (AKA parameter codes) to parameter names.
    *
    * Each parameter "code"/"index"/"ID number" is a positive integer that identifies the
    * given parameter in the input XML.
    *
    * The keys of this map need not form a contiguous sequence.
    * E.g. in the case where a parameter is deleted altogether from the simulation code.
    */
    const std::unordered_map<int, ParameterName> idCodeToNameMap
    {
        { 1,  ParameterName::NEG_LOG_ONE_MINUS_SINF },
        { 2,  ParameterName::E_STAR },
        { 3,  ParameterName::SIMM },
        { 4,  ParameterName::X_STAR_P },
        { 5,  ParameterName::GAMMA_P },
        { 6,  ParameterName::SIGMA_I_SQ },
        { 7,  ParameterName::CUMULATIVE_Y_STAR },
        { 8,  ParameterName::CUMULATIVE_H_STAR },
        { 9,  ParameterName::NEG_LOG_ONE_MINUS_ALPHA_M },
        { 10, ParameterName::DECAY_M },
        { 11, ParameterName::SIGMA0_SQ },
        { 12, ParameterName::X_NU_STAR },
        { 13, ParameterName::Y_STAR_SQ },
        { 14, ParameterName::ALPHA },
        { 15, ParameterName::DENSITY_BIAS_NON_GARKI },
        { 16, ParameterName::BASELINE_AVAILABILITY_SHAPE },
        { 17, ParameterName::LOG_ODDS_RATIO_CF_COMMUNITY },
        { 18, ParameterName::INDIRECT_RISK_COFACTOR },
        { 19, ParameterName::NON_MALARIA_INFANT_MORTALITY },
        { 20, ParameterName::DENSITY_BIAS_GARKI },
        { 21, ParameterName::SEVERE_MALARIA_THRESHHOLD },
        { 22, ParameterName::IMMUNITY_PENALTY },
        { 23, ParameterName::IMMUNE_EFFECTOR_DECAY },
        { 24, ParameterName::COMORBIDITY_INTERCEPT },
        { 25, ParameterName::Y_STAR_HALF_LIFE },
        { 26, ParameterName::Y_STAR_1 },
        { 27, ParameterName::ASEXUAL_IMMUNITY_DECAY },
        { 28, ParameterName::Y_STAR_0 },
        // 29 corresponds to a now-deprecated parameter.  If adding a new parameter,
        // don't use 29.  Since some old scenarios likely define a parameter value
        // with index 29 already.
        { 30, ParameterName::CRITICAL_AGE_FOR_COMORBIDITY },
        { 31, ParameterName::MUELLER_RATE_MULTIPLIER },
        { 32, ParameterName::MUELLER_DENSITY_EXPONENT },
        { 33, ParameterName::CFR_SCALE_FACTOR },
        { 34, ParameterName::MEAN_LOCAL_MAX_DENSITY },
        { 35, ParameterName::SD_LOCAL_MAX_DENSITY },
        { 36, ParameterName::MEAN_DIFF_POS_DAYS },
        { 37, ParameterName::SD_DIFF_POS_DAYS },
        { 38, ParameterName::CFR_NEG_LOG_ALPHA }
    };

    /*
    * Defines the map from parameter names to parameter (optional) values.
    *
    * The map values are made optional because many scenarios will specify
    * values for less than the full set of parameters.
    *
    * Any given ParameterName has a value (other than std::nullopt) iff
    * either the XML explicitly states a value for the parameter or the
    * XML explicitly states to use some model which itself contains the
    * given parameter among its preset values.
    */
    std::unordered_map<ParameterName, std::optional<double>> nameToValueMap;
};

}
#endif
