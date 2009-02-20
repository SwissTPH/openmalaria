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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

//
// This file contains all the names
//
#ifdef new
#undef new
#define _REDEF_NEW
#endif
#include <xercesc/dom/DOM.hpp>

#undef new
#ifdef _REDEF_NEW
#define new DEBUG_NEW
#undef _REDEF_NEW
#endif 

XERCES_CPP_NAMESPACE_USE

#ifndef NAMES_H
#define NAMES_H

#define s_MODE "mode"
#define s_ASSIM_MODE "assimMode"
#define s_RELEASE "release"
#define s_WU_ID "wuID"
#define s_PARAMETERS "parameters"
#define s_LOWER_BOUND "lowerbound"
#define s_UPPER_BOUND "upperbound"
#define s_POP_PERCENT "poppercent"
#define s_TIME "time"
#define s_MDA "MDA"
#define s_IRS "IRS"
#define s_VACCINATE "vaccinate"
#define s_CHANGE_EIR "changeEIR"
#define s_CHANGE_HS "changeHS"
#define s_ITN "ITN"
#define s_IPT "IPT"
#define s_TIMED "timed"
#define s_MIN_AGE "minAge"
#define s_MAX_AGE "maxAge"
#define s_COVERAGE "coverage"
#define s_SURVEYS "surveys"
#define s_AGE_GROUP "ageGroup"
#define s_ANALYSIS_NO "analysisNo"
#define s_NAME "name"
#define s_MAXIMUM_AGE_YEARS "maximumAgeYrs"
#define s_MINIMUM_AGE_YEARS "minimumAgeYrs"
#define s_POP_SIZE "popSize"
#define s_SIMULATION_DURATION "simulationDuration"
#define s_DEMOGRAPHY "demography" 
#define s_MONITORING "monitoring"
#define s_INTERVENTIONS "interventions"
#define s_ENTO_DATA "entoData"
#define s_DETECTION_LIMIT "detectionLimit"
#define s_SUMMARY_OPTION "summaryOption"

#define s_INPUTTYPE "inputType"
#define s_EIR "EIR"
#define s_ANOPHELES "anopheles"
#define s_USENV0GUESS "useNv0Guess"

#define s_MAX_AGE_YRS "maxAgeYrs"
#define s_MIN_AGE_YRS "minAgeYrs"
#define s_MODEL_VERSION "modelVersion"
#define s_HEALTH_SYSTEM "healthSystem"
#define s_HEALTH_SYSTEM_MEMORY "healthSystemMemory"
#define s_DRUG_REGIMEN "drugRegimen"
#define s_INITIAL_ACR "initialACR"
#define s_COMPLIANCE "compliance"
#define s_NON_COMPLIERS_EFFECTIVE "nonCompliersEffective"
#define s_P_SEEK_OFFICIAL_CARE_UNCOMPLICATED_1 "pSeekOfficialCareUncomplicated1"
#define s_P_SEEK_OFFICIAL_CARE_UNCOMPLICATED_2 "pSeekOfficialCareUncomplicated2"
#define s_P_SELF_TREAT_UNCOMPLICATED "pSelfTreatUncomplicated"
#define s_P_SEEK_OFFICIAL_CARE_SEVERE "pSeekOfficialCareSevere"
#define s_P_SEQUELAE_INPATIENT "pSequelaeInpatient"
#define s_FIRST_LINE "firstLine"
#define s_SECOND_LINE "secondLine"
#define s_INPATIENT "inpatient"
#define s_CQ "CQ"
#define s_SP "SP"
#define s_AQ "AQ"
#define s_SPAQ "SPAQ"
#define s_ACT "ACT"
#define s_QN "QN"
#define s_SELF_TREATMENT "selfTreatment"
#define s_VACCINE_DESCRIPTION "vaccineDescription"
#define s_VACCINE_TYPE "vaccineType"
#define s_HALF_LIFE_YRS "halfLifeYrs"
#define s_EFFICACY_B "efficacyB"
#define s_INITIAL_EFFICACY "initialEfficacy"
#define s_ITN_DESCRIPTION "ITNdescription"
#define s_IPT_DESCRIPTION "IPTdescription"
#define s_PU_0 "pu0"
#define s_PU_1 "pu1"
#define s_SPOROGONY_GONOTROPHY "sporogonyGonotrophy"
#define s_CONTINUOUS "continuous"
#define s_VACCINE "vaccine"
#define s_TARGET_AGE_YRS "targetAgeYrs"
#define s_BEST "best"
#define s_DESCRIPTION "description"
#define s_UNITS "units"
#define s_ASSUMPTION "assumption"
#define s_SOURCES "sources"
#define s_VALUE "value"
#define s_CFR "CFR"
#define s_CFRValue "cfr"
#define s_GROUP "group"
#define s_ITEM "item"
#define s_LATENT_P "latentp"
#define s_DELTA "delta"
#define s_NSPORE "nspore"
#define s_INTERVAL "interval"
#define s_I_SEED "iseed"
#define s_PARAMETER "parameter"
#define s_NUMBER "number" 
#define s_GROWTHRATE "growthrate"

#define s_IPTIDESCRIPTION "iptiDescription"
#define s_IPTIEFFECT "iptiEffect"
#define s_INFGENOTYPE "infGenotype"
#define s_FREQ "freq"
#define s_ACR "ACR"
#define s_PROPH "proph"
#define s_TOLPERIOD "tolPeriod"
#define s_ATTEN "atten"
#define s_IPTI "ipti"

#define s_PROTEOME "proteome"
#define s_CONTENT "content"
#define s_DISTRIBUTION "distribution"
#define s_POSITION "position"
#define s_CNV "cnv"
#define s_AMINOS "aminos"
#define s_PROPORTION "proportion"
#define s_FITNESS "fitness"

#define s_CASE_MANAGEMENT "caseManagement"
#define s_CASE_MANAGEMENTS "caseManagements"
#define s_DECISIONS "decisions"
#define s_DECISION "decision"
#define s_QTY "qty"
#define s_P "p"
#define s_ID "id"
#define s_MEDICATE "medicate"
#define s_UC1 "uc1"
#define s_UC2 "uc2"
#define s_SEV "sev"
#define s_NMF "nmf"
#define s_ENDPOINT "endPoint"

#endif
