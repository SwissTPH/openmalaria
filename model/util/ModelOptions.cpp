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

#include "util/ModelOptions.h"
#include "util/CommandLine.h"
#include "util/errors.h"
#include "inputData.h"

#include <sstream>
#include <iostream>

namespace OM { namespace util {
    uint32_t ModelOptions::optSet;
    
    // Utility: converts option strings to codes and back
    class OptionCodeMap {
	// Lookup table to translate the strings used in the XML file to the internal enumerated values:
	map<string,OptionCodes> codeMap;
	
    public:
	OptionCodeMap () {
	    codeMap["PENALISATION_EPISODES"] = PENALISATION_EPISODES;
	    codeMap["NEGATIVE_BINOMIAL_MASS_ACTION"] = NEGATIVE_BINOMIAL_MASS_ACTION;
	    codeMap["ATTENUATION_ASEXUAL_DENSITY"] = ATTENUATION_ASEXUAL_DENSITY;
	    codeMap["LOGNORMAL_MASS_ACTION"] = LOGNORMAL_MASS_ACTION;
	    codeMap["NO_PRE_ERYTHROCYTIC"] = NO_PRE_ERYTHROCYTIC;
	    codeMap["MAX_DENS_CORRECTION"] = MAX_DENS_CORRECTION;
	    codeMap["INNATE_MAX_DENS"] = INNATE_MAX_DENS;
	    // 	codeMap["MAX_DENS_RESET"] = MAX_DENS_RESET;
	    codeMap["DUMMY_WITHIN_HOST_MODEL"] = DUMMY_WITHIN_HOST_MODEL;
	    codeMap["PREDETERMINED_EPISODES"] = PREDETERMINED_EPISODES;
	    codeMap["NON_MALARIA_FEVERS"] = NON_MALARIA_FEVERS;
	    codeMap["INCLUDES_PK_PD"] = INCLUDES_PK_PD;
	    codeMap["CLINICAL_EVENT_SCHEDULER"] = CLINICAL_EVENT_SCHEDULER;
	    codeMap["MUELLER_PRESENTATION_MODEL"] = MUELLER_PRESENTATION_MODEL;
	    codeMap["TRANS_HET"] = TRANS_HET;
	    codeMap["COMORB_HET"] = COMORB_HET;
	    codeMap["TREAT_HET"] = TREAT_HET;
	    codeMap["COMORB_TRANS_HET"] = COMORB_TRANS_HET;
	    codeMap["TRANS_TREAT_HET"] = TRANS_TREAT_HET;
	    codeMap["COMORB_TREAT_HET"] = COMORB_TREAT_HET;
	    codeMap["TRIPLE_HET"] = TRIPLE_HET;
	    codeMap["EMPIRICAL_WITHIN_HOST_MODEL"] = EMPIRICAL_WITHIN_HOST_MODEL;
            codeMap["MOLINEAUX_WITHIN_HOST_MODEL"] = MOLINEAUX_WITHIN_HOST_MODEL;
            codeMap["PENNY_WITHIN_HOST_MODEL"] = PENNY_WITHIN_HOST_MODEL;
            codeMap["GARKI_DENSITY_BIAS"] = GARKI_DENSITY_BIAS;
            codeMap["IPTI_SP_MODEL"] = IPTI_SP_MODEL;
            codeMap["REPORT_ONLY_AT_RISK"] = REPORT_ONLY_AT_RISK;
	    
	    codeMap["MEAN_DURATION_GAMMA"] = MEAN_DURATION_GAMMA;
	    codeMap["FIRST_LOCAL_MAXIMUM_GAMMA"]=FIRST_LOCAL_MAXIMUM_GAMMA;
	    codeMap["IMMUNE_THRESHOLD_GAMMA"]=IMMUNE_THRESHOLD_GAMMA;
	    codeMap["UPDATE_DENSITY_GAMMA"]=UPDATE_DENSITY_GAMMA;
	    codeMap["PARASITE_REPLICATION_GAMMA"]=PARASITE_REPLICATION_GAMMA;
            codeMap["VECTOR_LIFE_CYCLE_MODEL"]=VECTOR_LIFE_CYCLE_MODEL;
	}
	
	OptionCodes operator[] (const string s) {
	    map<string,OptionCodes>::iterator codeIt = codeMap.find (s);
	    if (codeIt == codeMap.end()) {
		ostringstream msg;
		msg << "Unrecognised model option: ";
		msg << s;
		throw xml_scenario_error(msg.str());
	    }
	    return codeIt->second;
	}
	// reverse-lookup in map; only used for error/debug printing so efficiency is unimportant
	// doesn't ensure code is unique in the map either
	string toString (const OptionCodes code) {
	    for (map<string,OptionCodes>::iterator codeIt = codeMap.begin(); codeIt != codeMap.end(); ++codeIt) {
		if (codeIt->second == code)
		    return codeIt->first;
	    }
	    throw TRACED_EXCEPTION_DEFAULT ("toString called with unknown code");	// this is a code error
	}
    };
    
    void ModelOptions::init () {
	OptionCodeMap codeMap;
	
	// State of all default options:
	bitset<NUM_OPTIONS> defaultOptSet;
	defaultOptSet.set (MAX_DENS_CORRECTION);
	
	// Set optSet to defaults, then override any given in the XML file:
	bitset<NUM_OPTIONS> optSet_bs = defaultOptSet;
	
	const scnXml::OptionSet::OptionSequence& optSeq = InputData().getModel().getModelOptions().getOption();
	for (scnXml::OptionSet::OptionConstIterator it = optSeq.begin(); it != optSeq.end(); ++it) {
	    optSet_bs[codeMap[it->getName()]] = it->getValue();
	}
	
	// Print non-default model options:
	if (CommandLine::option (CommandLine::PRINT_MODEL_OPTIONS)) {
	    cout << "Non-default model options:";
	    for (int i = 0; i < NUM_OPTIONS; ++i) {
		if (optSet_bs[i] != defaultOptSet[i])
		    cout << "\t" << codeMap.toString(OptionCodes(i)) << "=" << optSet_bs[i];
	    }
	    cout << endl;
	}
	
	// Test for incompatible options
	
	// An incompatibility triangle, listing options incompatible with each option
	// Doesn't check required versions
	// Originally from "description of variables for interface" excel sheet
	bitset<NUM_OPTIONS> incompatibilities[NUM_OPTIONS];	// all default to 0
	
	incompatibilities[NEGATIVE_BINOMIAL_MASS_ACTION]
	    .set(LOGNORMAL_MASS_ACTION)
	    .set(TRANS_HET)	.set(COMORB_TRANS_HET)
	    .set(TRANS_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[LOGNORMAL_MASS_ACTION]
	    .set(TRANS_HET)	.set(COMORB_TRANS_HET)
	    .set(TRANS_TREAT_HET)	.set(TRIPLE_HET);
	
	// Note: MAX_DENS_CORRECTION is irrelevant when using new
	// within-host models, but we don't mark it incompatible so that we can
	// leave MAX_DENS_CORRECTION on by default.
	incompatibilities[DUMMY_WITHIN_HOST_MODEL]
            .set(PENNY_WITHIN_HOST_MODEL)
	    .set(EMPIRICAL_WITHIN_HOST_MODEL)
	    .set(MOLINEAUX_WITHIN_HOST_MODEL)
	    .set(MEAN_DURATION_GAMMA)
	    .set(FIRST_LOCAL_MAXIMUM_GAMMA)
	    .set(PARASITE_REPLICATION_GAMMA)
	    .set(IMMUNE_THRESHOLD_GAMMA)
	    .set(UPDATE_DENSITY_GAMMA);

	incompatibilities[EMPIRICAL_WITHIN_HOST_MODEL]
            .set(MOLINEAUX_WITHIN_HOST_MODEL)
            .set(PENNY_WITHIN_HOST_MODEL)
	    .set(MEAN_DURATION_GAMMA)
	    .set(FIRST_LOCAL_MAXIMUM_GAMMA)
	    .set(PARASITE_REPLICATION_GAMMA)
	    .set(IMMUNE_THRESHOLD_GAMMA)
	    .set(UPDATE_DENSITY_GAMMA);
	    
	    
	incompatibilities[MOLINEAUX_WITHIN_HOST_MODEL]
            .set(PENNY_WITHIN_HOST_MODEL)
	    .set(IMMUNE_THRESHOLD_GAMMA)
	    .set(UPDATE_DENSITY_GAMMA);
	    
	incompatibilities[PENNY_WITHIN_HOST_MODEL]
	    .set(MEAN_DURATION_GAMMA)
	    .set(FIRST_LOCAL_MAXIMUM_GAMMA)
	    .set(PARASITE_REPLICATION_GAMMA);
	    
	incompatibilities[NON_MALARIA_FEVERS]
	    .set(MUELLER_PRESENTATION_MODEL);
	
	incompatibilities[TRANS_HET]
	    .set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_HET]
	    .set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[TREAT_HET]
	    .set(COMORB_TRANS_HET)	.set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_TRANS_HET]
	    .set(TRANS_TREAT_HET)
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[TRANS_TREAT_HET]
	    .set(COMORB_TREAT_HET)	.set(TRIPLE_HET);
	incompatibilities[COMORB_TREAT_HET]
	    .set(TRIPLE_HET);
        
	for (size_t i = 0; i < NUM_OPTIONS; ++i) {
	    if (optSet_bs [i] && (optSet_bs & incompatibilities[i]).any()) {
		ostringstream msg;
		msg << "Incompatible model options: " << codeMap.toString(OptionCodes(i)) << "=" << optSet_bs[i]
			<< " is incompatible with flags:";
		bitset<NUM_OPTIONS> incompat = (optSet_bs & incompatibilities[i]);
		for (int j = 0; j < NUM_OPTIONS; ++j) {
		    if (incompat[j])
			msg << "\t" << codeMap.toString(OptionCodes(j)) << "=" << optSet_bs[j];
		}
		throw xml_scenario_error (msg.str());
	    }
	}
	
        // Required options (above table can't check these):
        if (optSet_bs[INNATE_MAX_DENS] && !optSet_bs[MAX_DENS_CORRECTION])
            throw xml_scenario_error ("INNATE_MAX_DENS requires MAX_DENS_CORRECTION");
        
        if( TimeStep::interval != 1 ){
            bitset<NUM_OPTIONS> require1DayTS;
            require1DayTS
                .set( DUMMY_WITHIN_HOST_MODEL )
                .set( INCLUDES_PK_PD )
                .set( CLINICAL_EVENT_SCHEDULER )
                .set( EMPIRICAL_WITHIN_HOST_MODEL )
                .set( MOLINEAUX_WITHIN_HOST_MODEL )
                .set( PENNY_WITHIN_HOST_MODEL );
            
            for (size_t i = 0; i < NUM_OPTIONS; ++i) {
                if (optSet_bs [i] && require1DayTS[i]) {
                    ostringstream msg;
                    msg << "Model option " << codeMap.toString(OptionCodes(i)) << " is only compatible with a 1-day timestep.";
                    throw xml_scenario_error (msg.str());
                }
            }
        }
        if( TimeStep::interval != 5){
            bitset<NUM_OPTIONS> require5DayTS;
            require5DayTS
                .set( ATTENUATION_ASEXUAL_DENSITY )
                /* This is only relevant to 5-day, but enabled by default.
                .set( MAX_DENS_CORRECTION ) */
                .set( INNATE_MAX_DENS )
                .set( IPTI_SP_MODEL )
                .set( REPORT_ONLY_AT_RISK );
            
            for (size_t i = 0; i < NUM_OPTIONS; ++i) {
                if (optSet_bs [i] && require5DayTS[i]) {
                    ostringstream msg;
                    msg << "Model option " << codeMap.toString(OptionCodes(i)) << " is only compatible with a 5-day timestep.";
                    throw xml_scenario_error (msg.str());
                }
            }
        }
	
	// Convert from bitset to more performant integer with binary operations
	// Note: use bitset up to now to restrict use of binary operators to
	// where it has significant performance benefits.
	optSet = 0;
	for (size_t i = 0; i < NUM_OPTIONS; ++i) {
	    if (optSet_bs[i])
		optSet |= (1<<i);
	}
    }
    
} }
