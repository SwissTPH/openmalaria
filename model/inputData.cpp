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


#include "inputData.h"
#include "util/BoincWrapper.h"
#include "util/errors.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <boost/format.hpp>

namespace OM {
    using boost::format;
    
// Initialization functions:

void InputDataType::initParameterValues()
{
    // set parameters
    const scnXml::Parameters::ParameterSequence& paramSeq = scenario->getModel().getParameters().getParameter();
    for (scnXml::Parameters::ParameterConstIterator it = paramSeq.begin(); it != paramSeq.end(); ++it) {
        int i = it->getNumber();
        if (i < 0 || i >= Params::MAX)
	    throw util::xml_scenario_error( (format("parameter with invalid index %1%") %i).str() );
	if( !parameterValues.insert( make_pair( i, it->getValue() ) ).second )
	    throw util::xml_scenario_error( (format("parameter with index %1% described twice") %i).str() );
    }
}

void InputDataType::initTimedInterventions()
{
    if (scenario->getInterventions().getContinuous().present()) {
        const scnXml::Continuous& contI = scenario->getInterventions().getContinuous().get();
        if (contI.getVaccine().size())
            activeInterventions.set (Interventions::VACCINE, true);
	if (contI.getITN().size())
	    activeInterventions.set (Interventions::ITN, true);
	if (contI.getIpti().size())
	    activeInterventions.set (Interventions::IPTI, true);
	if (contI.getCohort().size())
	    activeInterventions.set (Interventions::COHORT, true);
    }

    if (scenario->getInterventions().getTimed().present()) {
        const scnXml::Timed::InterventionSequence& interventionSeq =
            scenario->getInterventions().getTimed().get().getIntervention();
        for (scnXml::Timed::InterventionConstIterator it (interventionSeq.begin()); it != interventionSeq.end(); ++it) {
            TimeStep time = TimeStep(it->getTime());
            if (timedInterventions.count (time)) {
                ostringstream msg;
                msg << "Error: multiple timed interventions with time: " << time;
                throw util::xml_scenario_error (msg.str());
            }
            timedInterventions[time] = & (*it);

            if (it->getChangeHS().present())
                activeInterventions.set (Interventions::CHANGE_HS, true);
            if (it->getChangeEIR().present())
                activeInterventions.set (Interventions::CHANGE_EIR, true);
            if (it->getVaccinate().present())
                activeInterventions.set (Interventions::VACCINE, true);
            if (it->getMDA().present())
                activeInterventions.set (Interventions::MDA, true);
            if (it->getIpti().present())
                activeInterventions.set (Interventions::IPTI, true);
            if (it->getITN().present())
                activeInterventions.set (Interventions::ITN, true);
            if (it->getIRS().present())
                activeInterventions.set (Interventions::IRS, true);
            if (it->getVectorAvailability().present())
                activeInterventions.set (Interventions::VEC_AVAIL, true);
            if (it->getImmuneSuppression().present())
                activeInterventions.set (Interventions::IMMUNE_SUPPRESSION, true);
            if (it->getCohort().present())
                activeInterventions.set (Interventions::COHORT, true);
            if (it->getLarviciding().present())
                activeInterventions.set (Interventions::LARVICIDING, true);
            if (it->getInsertR_0Case().present()){
                activeInterventions.set (Interventions::R_0_CASE, true);
		// uses vaccines but see note in Vaccine::initParameters()
                // activeInterventions.set (Interventions::VACCINE, true);
	    }
            if (it->getImportedInfectionsPerThousandHosts().present())
            	activeInterventions.set (Interventions::IMPORTED_INFECTIONS, true);
            if (it->getUninfectVectors().present())
                activeInterventions.set (Interventions::UNINFECT_VECTORS, true);

        }
    }
}


util::Checksum InputDataType::createDocument (std::string lXmlFile)
{
    xmlFileName = lXmlFile;
    //Parses the document
    
    // Opening by filename causes a schema lookup in the scenario file's dir,
    // which no longer works (schemas moved).
    // Opening with a stream causes it to look in the working directory.
    //NOTE: it'd be nice if this used Global::lookupResource for the schema.
    
    // Note that the schema location can be set manually by passing properties,
    // but we won't necessarily have the right schema version associated with
    // the XML file in that case.
    ifstream fileStream (lXmlFile.c_str(), ios::binary);
    if (!fileStream.good()){
	string msg = "Error: unable to open "+lXmlFile;
	throw util::xml_scenario_error (msg);
    }
    scenario = scnXml::parseScenario (fileStream);
    util::Checksum cksum = util::Checksum::generate (fileStream);
    fileStream.close ();
    int scenarioVersion = scenario->getSchemaVersion();
    if (scenarioVersion < SCHEMA_VERSION) {
	ostringstream msg;
	msg<<lXmlFile<<" uses "
	    << ((scenarioVersion < SCHEMA_VERSION_OLDEST_COMPATIBLE) ? "an" : "a potentially")
	    <<" incompatible old schema version ("<<scenarioVersion<<"; current is "
	    <<SCHEMA_VERSION<<"). Use SchemaTranslator to update.";
	if (scenarioVersion < SCHEMA_VERSION_OLDEST_COMPATIBLE) {
	    throw util::xml_scenario_error (msg.str());
	} else {
	    cerr<<"Warning: "<<msg.str()<<endl;
	}
    }
    if (scenarioVersion > SCHEMA_VERSION)
        throw util::xml_scenario_error ("Error: new schema version unsupported");

    initParameterValues();
    initTimedInterventions();
    return cksum;
}

void InputDataType::saveDocument()
{
    if (documentChanged) {
        // get the "basename" (file name without path) of xmlFileName as a C string:
        const char* lastFS = strrchr (xmlFileName.c_str(), '/');
        const char* lastBS = strrchr (xmlFileName.c_str(), '\\');
        const char* baseName = lastBS > lastFS ? lastBS : lastFS;
        if (baseName == NULL) // no path separator found; use whole string
            baseName = xmlFileName.c_str();
        else
            ++baseName;  // start at next character
	
        ofstream outStream (baseName);
	// Set schema file. Unfortunately we don't know what it was in input
	// file, so this is only a guess.
        ostringstream schema;
        schema << "scenario_" << SCHEMA_VERSION << ".xsd";

        xml_schema::NamespaceInfomap map;
        map[""].name = "";
        map[""].schema = schema.str();
        scnXml::serializeScenario (outStream, *scenario, map);

        outStream.close();
    }
}

void InputDataType::freeDocument(){
    // Destructors should handle cleanup
}


double InputDataType::getParameter (size_t i)
{
    std::map<int, double>::const_iterator it = parameterValues.find( i );
    if( it == parameterValues.end() )
	throw util::xml_scenario_error( (format("parameter %1% required but not described") %i).str() );
    return it->second;
}
const scnXml::Intervention* InputDataType::getInterventionByTime (TimeStep time)
{
    std::map<TimeStep, const scnXml::Intervention*>::iterator i = timedInterventions.find (time);
    if (i != timedInterventions.end())
        return i->second;
    else
        return NULL;
}



InputDataType InputData;
}
